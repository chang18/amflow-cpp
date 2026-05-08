// SPDX-License-Identifier: MIT
// ibp::kira_run + kira_read_masters + kira_read_target_table.
//

#include "amflow/ibp/kira.hpp"

#include <cctype>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

extern "C" char** environ;

namespace amflow::ibp {

namespace fs = std::filesystem;

double kira_run(const KiraConfig& cfg,
                const std::string& dir,
                const std::string& log_path) {
    if (cfg.kira_executable.empty()) {
        throw std::invalid_argument("kira_run: kira_executable not set");
    }
    if (!fs::exists(cfg.kira_executable)) {
        throw std::runtime_error(
            "kira_run: kira executable not found at " + cfg.kira_executable);
    }
    if (!cfg.fermat_executable.empty() && !fs::exists(cfg.fermat_executable)) {
        throw std::runtime_error(
            "kira_run: fermat executable not found at " + cfg.fermat_executable);
    }
    if (!fs::is_directory(dir)) {
        throw std::runtime_error("kira_run: dir not found " + dir);
    }

    std::vector<std::string> args;
    args.push_back(cfg.kira_executable);
    args.push_back("-p" + std::to_string(cfg.n_thread));
    args.push_back("jobs.yaml");
    for (const auto& [name, value] : cfg.numeric_values) {
        if (name == "eps") {
            args.push_back("-sd=4-2*" + value);
        } else {
            args.push_back("-s" + name + "=" + value);
        }
    }

    std::vector<char*> argv_c;
    argv_c.reserve(args.size() + 1);
    for (auto& s : args) argv_c.push_back(const_cast<char*>(s.c_str()));
    argv_c.push_back(nullptr);

    std::vector<std::string> env_strings;
    bool fermat_set = false;
    for (char** e = environ; *e; ++e) {
        std::string entry(*e);
        if (entry.rfind("FERMATPATH=", 0) == 0) {
            if (!cfg.fermat_executable.empty()) {
                entry = "FERMATPATH=" + cfg.fermat_executable;
                fermat_set = true;
            }
        }
        env_strings.push_back(entry);
    }
    if (!fermat_set && !cfg.fermat_executable.empty()) {
        env_strings.push_back("FERMATPATH=" + cfg.fermat_executable);
    }
    std::vector<char*> envp_c;
    envp_c.reserve(env_strings.size() + 1);
    for (auto& s : env_strings) envp_c.push_back(const_cast<char*>(s.c_str()));
    envp_c.push_back(nullptr);

    int log_fd = -1;
    if (!log_path.empty()) {
        log_fd = ::open(log_path.c_str(),
                         O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (log_fd < 0) {
            throw std::runtime_error(
                "kira_run: cannot open log file " + log_path);
        }
    }

    auto t0 = std::chrono::steady_clock::now();

    pid_t pid = fork();
    if (pid < 0) {
        if (log_fd >= 0) close(log_fd);
        throw std::runtime_error(
            "kira_run: fork failed (" + std::string(strerror(errno)) + ")");
    }
    if (pid == 0) {
        if (chdir(dir.c_str()) != 0) {
            _exit(127);
        }
        if (log_fd >= 0) {
            dup2(log_fd, STDOUT_FILENO);
            dup2(log_fd, STDERR_FILENO);
            close(log_fd);
        }
        execve(cfg.kira_executable.c_str(), argv_c.data(), envp_c.data());
        _exit(127);
    }
    if (log_fd >= 0) close(log_fd);

    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        throw std::runtime_error(
            "kira_run: waitpid failed (" + std::string(strerror(errno)) + ")");
    }

    auto t1 = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration<double>(t1 - t0).count();

    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        throw std::runtime_error(
            "kira_run: kira exited with non-zero status (rc=" +
            std::to_string(WEXITSTATUS(status)) + ")");
    }
    return elapsed;
}

namespace {

void skip_ws(const std::string& s, std::size_t& i) {
    while (i < s.size() && std::isspace((unsigned char)s[i])) ++i;
}

std::string parse_ident(const std::string& s, std::size_t& i) {
    std::string out;
    while (i < s.size() &&
            (std::isalnum((unsigned char)s[i]) || s[i] == '_')) {
        out += s[i];
        ++i;
    }
    return out;
}

std::optional<qft::JIntegral>
parse_jintegral_inplace(const std::string& s, std::size_t& i) {
    skip_ws(s, i);
    std::string fam = parse_ident(s, i);
    if (fam.empty()) return std::nullopt;
    skip_ws(s, i);
    if (i >= s.size() || (s[i] != '(' && s[i] != '[')) return std::nullopt;
    char closer = (s[i] == '(') ? ')' : ']';
    ++i;
    std::vector<long> indices;
    while (true) {
        skip_ws(s, i);
        long sign = 1;
        if (i < s.size() && s[i] == '-') { sign = -1; ++i; }
        else if (i < s.size() && s[i] == '+') { ++i; }
        std::string num;
        while (i < s.size() && std::isdigit((unsigned char)s[i])) {
            num += s[i]; ++i;
        }
        if (num.empty()) return std::nullopt;
        indices.push_back(sign * std::stol(num));
        skip_ws(s, i);
        if (i < s.size() && s[i] == ',') {
            ++i; continue;
        }
        if (i < s.size() && s[i] == closer) {
            ++i; break;
        }
        return std::nullopt;
    }
    return qft::JIntegral(fam, indices);
}

}  // namespace

std::vector<qft::JIntegral>
kira_read_masters(const qft::FamilyConfig& fc, const std::string& dir) {
    std::string path = dir + "/results/" + fc.family + "/masters";
    if (!fs::exists(path)) {
        return {};
    }
    std::ifstream f(path);
    if (!f) throw std::runtime_error("kira_read_masters: cannot open " + path);

    std::vector<qft::JIntegral> out;
    std::string line;
    while (std::getline(f, line)) {
        std::size_t i = 0;
        skip_ws(line, i);
        if (i >= line.size()) continue;
        if (line[i] == '#') continue;
        auto j = parse_jintegral_inplace(line, i);
        if (!j.has_value()) {
            throw std::runtime_error(
                "kira_read_masters: cannot parse line: " + line);
        }
        out.push_back(std::move(*j));
    }
    return out;
}

namespace {

std::string slurp(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("slurp: cannot open " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::size_t find_matching_close(const std::string& s, std::size_t i,
                                  char open, char close) {
    long depth = 0;
    for (; i < s.size(); ++i) {
        if (s[i] == open) ++depth;
        else if (s[i] == close) {
            if (--depth == 0) return i;
        }
    }
    throw std::runtime_error("find_matching_close: unmatched");
}

std::vector<std::string>
split_top_level(const std::string& s, char sep) {
    std::vector<std::string> out;
    long depth = 0;
    std::size_t start = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '(' || c == '{' || c == '[') ++depth;
        else if (c == ')' || c == '}' || c == ']') --depth;
        else if (depth == 0 && c == sep) {
            out.push_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    out.push_back(s.substr(start));
    return out;
}

std::string trim(const std::string& s) {
    std::size_t a = 0;
    while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
    std::size_t b = s.size();
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

std::size_t find_top_arrow(const std::string& s) {
    long depth = 0;
    for (std::size_t i = 0; i + 1 < s.size(); ++i) {
        char c = s[i];
        if (c == '(' || c == '{' || c == '[') ++depth;
        else if (c == ')' || c == '}' || c == ']') --depth;
        else if (depth == 0 && c == '-' && s[i + 1] == '>') {
            return i;
        }
    }
    return std::string::npos;
}

std::vector<std::string> split_terms(const std::string& s) {
    std::vector<std::string> out;
    long depth = 0;
    std::size_t start = 0;
    bool first = true;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '(' || c == '{' || c == '[') ++depth;
        else if (c == ')' || c == '}' || c == ']') --depth;
        else if (depth == 0 && (c == '+' || c == '-')) {
            if (i == 0 || first) {
                first = false;
                continue;
            }
            std::size_t k = i;
            while (k > 0 && std::isspace((unsigned char)s[k - 1])) --k;
            if (k == 0) continue;
            char prev = s[k - 1];
            if (prev == '*' || prev == '/' || prev == '^' || prev == '(' ||
                prev == ',' || prev == '=' || prev == '+' || prev == '-') {
                continue;
            }
            out.push_back(s.substr(start, i - start));
            start = i;
            first = false;
        } else if (!std::isspace((unsigned char)c)) {
            first = false;
        }
    }
    out.push_back(s.substr(start));
    return out;
}

struct TermParts {
    std::string coef;
    std::optional<qft::JIntegral> j;
};

TermParts split_term_into_coef_and_jintegral(
    const std::string& term, const std::string& family) {
    TermParts r;
    std::size_t pos = std::string::npos;
    char open_char = '(';
    char close_char = ')';
    {
        std::size_t p1 = term.find(family + "(");
        std::size_t p2 = term.find(family + "[");
        std::size_t cand;
        if (p1 == std::string::npos && p2 == std::string::npos) {
            r.coef = trim(term);
            return r;
        }
        if (p1 == std::string::npos) cand = p2;
        else if (p2 == std::string::npos) cand = p1;
        else cand = std::min(p1, p2);
        pos = cand;
        if (pos < term.size() && term[pos + family.size()] == '[') {
            open_char = '['; close_char = ']';
        }
        if (pos > 0) {
            char p = term[pos - 1];
            if (std::isalnum((unsigned char)p) || p == '_') {
                r.coef = trim(term);
                return r;
            }
        }
    }
    std::size_t open = pos + family.size();
    std::size_t close = find_matching_close(term, open, open_char, close_char);
    std::size_t i = pos;
    auto j = parse_jintegral_inplace(term, i);
    if (!j.has_value() || i != close + 1) {
        r.coef = trim(term);
        return r;
    }
    r.j = std::move(*j);

    auto strip_mul_space = [](std::string& s, bool from_back) {
        if (from_back) {
            while (!s.empty() && (s.back() == '*' ||
                    std::isspace((unsigned char)s.back()))) s.pop_back();
        } else {
            while (!s.empty() && (s.front() == '*' ||
                    std::isspace((unsigned char)s.front())))
                s.erase(s.begin());
        }
    };
    std::string prefix = trim(term.substr(0, pos));
    std::string suffix = trim(term.substr(close + 1));
    strip_mul_space(prefix, true);
    strip_mul_space(suffix, false);
    if (!prefix.empty() && prefix.front() == '+') {
        prefix = trim(prefix.substr(1));
    }

    std::string coef;
    if (prefix.empty() && suffix.empty()) {
        coef = "1";
    } else if (prefix.empty()) {
        coef = suffix;
    } else if (suffix.empty()) {
        coef = prefix;
    } else {
        coef = "(" + prefix + ")*(" + suffix + ")";
    }
    if (coef == "-") coef = "-1";
    r.coef = std::move(coef);
    return r;
}

}  // namespace

std::vector<KiraReductionRule>
kira_read_target_table(const qft::FamilyConfig& fc, const std::string& dir) {
    std::string path = dir + "/results/" + fc.family + "/kira_target.m";
    if (!fs::exists(path)) return {};
    std::string content = slurp(path);

    std::size_t open = content.find('{');
    if (open == std::string::npos) return {};
    std::size_t close = find_matching_close(content, open, '{', '}');
    std::string body = content.substr(open + 1, close - open - 1);

    auto rules_strs = split_top_level(body, ',');

    std::vector<KiraReductionRule> out;
    for (auto& rs : rules_strs) {
        std::string r = trim(rs);
        if (r.empty()) continue;
        std::size_t arrow = find_top_arrow(r);
        if (arrow == std::string::npos) continue;
        std::string lhs_str = trim(r.substr(0, arrow));
        std::string rhs_str = trim(r.substr(arrow + 2));

        std::size_t k = 0;
        auto lhs = parse_jintegral_inplace(lhs_str, k);
        if (!lhs.has_value()) {
            throw std::runtime_error(
                "kira_read_target_table: cannot parse LHS '" + lhs_str + "'");
        }

        KiraReductionRule kr;
        kr.lhs = std::move(*lhs);

        auto terms = split_terms(rhs_str);
        for (auto& t : terms) {
            std::string ts = trim(t);
            if (ts.empty()) continue;
            std::string ts_no_sign = ts;
            if (!ts_no_sign.empty() && (ts_no_sign[0] == '+' || ts_no_sign[0] == '-')) {
                ts_no_sign = trim(ts_no_sign.substr(1));
            }
            if (ts_no_sign == "0") continue;
            TermParts tp = split_term_into_coef_and_jintegral(ts, fc.family);
            if (!tp.j.has_value()) {
                throw std::runtime_error(
                    "kira_read_target_table: term has no J integral: '"
                    + ts + "'");
            }
            kr.rhs.emplace_back(std::move(tp.coef), std::move(*tp.j));
        }
        out.push_back(std::move(kr));
    }
    return out;
}

}  // namespace amflow::ibp
