(* 3-loop massless cut banana (4 cut propagators between two external points).
   Family: cutbanana[a, b, c, d, e, f, g, h, i] with 4 actual propagators
       l1^2, l2^2, l3^2, (l1+l2+l3-p)^2
   plus 5 ISPs (l1+p)^2, (l2+p)^2, (l3+p)^2, (l1-l2)^2, (l1-l3)^2 to
   complete the 9-SP basis (3 loops + 1 leg gives 9 independent SPs).
   Top sector j[cutbanana, 1, 1, 1, 1, 0, 0, 0, 0, 0] satisfies
   phase_volume_q (|var| = 4 = loops + 1, all live cuts == 1), so
   AMFSystemSetupMaster picks the Cutkosky branch at L = 3: cut is
   cleared, eta is injected on l1^2, and the final answer is
   Im[uncut]*prefactor with prefactor =
   2 (Pi^(2-eps) * (2 Pi)^(2 eps - 4))^L * (-1)^(L+1) at L = 3.
   We take s -> 4 (above the massless 4-particle threshold). epsorder
   pushed deep so the Laurent at eps=1/100 sums tightly enough to
   compare with C++ at ~1e-30 relative precision. *)
name        = "cutbanana_3L";
family      = cutbanana;
loops       = {l1, l2, l3};
legs        = {p};
conservation = {};
replacement = {p^2 -> s};
propagators = {l1^2, l2^2, l3^2, (l1 + l2 + l3 - p)^2,
               (l1 + p)^2, (l2 + p)^2, (l3 + p)^2,
               (l1 - l2)^2, (l1 - l3)^2};
cut         = {1, 1, 1, 1, 0, 0, 0, 0, 0};
numeric     = {s -> 4};
targets     = {j[cutbanana, 1, 1, 1, 1, 0, 0, 0, 0, 0]};
precision   = 50;
epsorder    = 22;
