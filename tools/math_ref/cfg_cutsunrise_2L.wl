(* 2-loop massless cut sunrise.
   Family: cutsunrise[a, b, c, d, e] with 3 actual propagators
       l1^2, l2^2, (l1+l2-p)^2
   plus 2 ISPs (l1+p)^2 and (l2+p)^2 to complete the 5-SP basis.
   Top sector j[cutsunrise, 1, 1, 1, 0, 0] satisfies phase_volume_q
   (|var| = 3 = loops + 1, all live cuts == 1), so AMFSystemSetupMaster
   picks the Cutkosky branch at L=2: cut is cleared, eta is injected,
   and the final answer is Im[uncut]*prefactor with prefactor =
   2 (Pi^(2-eps) * (2 Pi)^(2 eps - 4))^L * (-1)^(L+1) at L = 2.
   We take s -> 4 (above the massless 3-particle threshold). epsorder
   pushed to 14 so the Laurent expansion is accurate enough at
   eps = 1/100 to compare against C++ at ~1e-30 relative precision. *)
name        = "cutsunrise_2L";
family      = cutsunrise;
loops       = {l1, l2};
legs        = {p};
conservation = {};
replacement = {p^2 -> s};
propagators = {l1^2, l2^2, (l1 + l2 - p)^2, (l1 + p)^2, (l2 + p)^2};
cut         = {1, 1, 1, 0, 0};
numeric     = {s -> 4};
targets     = {j[cutsunrise, 1, 1, 1, 0, 0]};
precision   = 50;
epsorder    = 22;
