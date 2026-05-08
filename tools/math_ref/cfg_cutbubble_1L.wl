(* 1-loop massless cut bubble.
   Family: cutbubble[a, b] with propagators l^2 and (l+p)^2, both on cut.
   Top sector j[cutbubble, 1, 1] satisfies phase_volume_q (|var| = loops + 1
   with all cuts == 1), so AMFSystemSetupMaster picks the Cutkosky branch:
   the cut is cleared, eta is injected on a propagator, and the final
   answer is taken as Im[uncut]*prefactor with prefactor =
   2 (Pi^(2-eps) * (2 Pi)^(2 eps - 4))^L * (-1)^(L+1) at L = 1.
   We take s -> 4 (above threshold so the cut is real) and push epsorder
   to 6, plenty for eps=1/100 sampling. *)
name        = "cutbubble_1L";
family      = cutbubble;
loops       = {l};
legs        = {p};
conservation = {};
replacement = {p^2 -> s};
propagators = {l^2, (l + p)^2};
cut         = {1, 1};
numeric     = {s -> 4};
targets     = {j[cutbubble, 1, 1]};
precision   = 50;
epsorder    = 14;
