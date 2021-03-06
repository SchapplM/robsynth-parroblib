% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P6RRPRRR14V3G8A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [24x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P6RRPRRR14V3G8A0_convert_par2_MPV_fixb.m

% Output:
% taugX [6x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 23:36
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G8A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(24,1)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: qJ has to be [3x6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: Koppelpunkt has to be [6x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'P6RRPRRR14V3G8A0_gravload_para_pf_mdp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 1
% StartTime: 2020-03-12 23:36:10
% EndTime: 2020-03-12 23:36:12
% DurationCPUTime: 1.67s
% Computational Cost: add. (11151->1066), mult. (26343->1946), div. (936->12), fcn. (28167->66), ass. (0->732)
unknown=NaN(6,1);
t1 = 0.1e1 / qJ(3,1);
t2 = sin(qJ(1,1));
t3 = cos(legFrame(1,3));
t5 = cos(qJ(1,1));
t6 = sin(legFrame(1,3));
t8 = t3 * t2 + t6 * t5;
t9 = t8 * t1;
t10 = cos(legFrame(1,2));
t11 = sin(qJ(2,1));
t12 = 0.1e1 / t11;
t13 = t12 * t10;
t14 = cos(legFrame(1,1));
t16 = sin(legFrame(1,1));
t18 = -t14 * g(2) - g(3) * t16;
t20 = sin(legFrame(1,2));
t21 = g(2) * t20;
t23 = g(3) * t20;
t25 = t10 * g(1);
t26 = -t14 * t23 + t16 * t21 + t25;
t28 = t5 * t18 + t26 * t2;
t31 = -t18 * t2;
t32 = t5 * t26 + t31;
t34 = t3 * t28 + t6 * t32;
t37 = 0.1e1 / qJ(3,2);
t38 = sin(qJ(1,2));
t39 = cos(legFrame(2,3));
t41 = cos(qJ(1,2));
t42 = sin(legFrame(2,3));
t44 = t39 * t38 + t42 * t41;
t45 = t44 * t37;
t46 = cos(legFrame(2,2));
t47 = sin(qJ(2,2));
t48 = 0.1e1 / t47;
t49 = t48 * t46;
t50 = cos(legFrame(2,1));
t52 = sin(legFrame(2,1));
t54 = -t50 * g(2) - g(3) * t52;
t56 = sin(legFrame(2,2));
t57 = g(2) * t56;
t59 = t50 * t56;
t61 = t46 * g(1);
t62 = -g(3) * t59 + t52 * t57 + t61;
t64 = t62 * t38 + t41 * t54;
t67 = -t54 * t38;
t68 = t41 * t62 + t67;
t70 = t39 * t64 + t42 * t68;
t73 = 0.1e1 / qJ(3,3);
t74 = sin(qJ(1,3));
t75 = cos(legFrame(3,3));
t77 = cos(qJ(1,3));
t78 = sin(legFrame(3,3));
t80 = t75 * t74 + t78 * t77;
t81 = t80 * t73;
t82 = cos(legFrame(3,2));
t83 = sin(qJ(2,3));
t84 = 0.1e1 / t83;
t85 = t84 * t82;
t86 = cos(legFrame(3,1));
t88 = sin(legFrame(3,1));
t90 = -t86 * g(2) - g(3) * t88;
t92 = sin(legFrame(3,2));
t93 = g(2) * t92;
t95 = g(3) * t92;
t97 = t82 * g(1);
t98 = -t86 * t95 + t88 * t93 + t97;
t100 = t98 * t74 + t77 * t90;
t103 = -t90 * t74;
t104 = t77 * t98 + t103;
t106 = t75 * t100 + t78 * t104;
t109 = 0.1e1 / qJ(3,4);
t110 = sin(qJ(1,4));
t111 = cos(legFrame(4,3));
t113 = cos(qJ(1,4));
t114 = sin(legFrame(4,3));
t116 = t111 * t110 + t114 * t113;
t117 = t116 * t109;
t118 = cos(legFrame(4,2));
t119 = sin(qJ(2,4));
t120 = 0.1e1 / t119;
t121 = t120 * t118;
t122 = cos(legFrame(4,1));
t124 = sin(legFrame(4,1));
t126 = -t122 * g(2) - g(3) * t124;
t128 = sin(legFrame(4,2));
t134 = t124 * g(2) * t128 - g(3) * t122 * t128 + t118 * g(1);
t136 = t134 * t110 + t113 * t126;
t139 = -t126 * t110;
t140 = t113 * t134 + t139;
t142 = t111 * t136 + t114 * t140;
t145 = 0.1e1 / qJ(3,5);
t146 = sin(qJ(1,5));
t147 = cos(legFrame(5,3));
t149 = cos(qJ(1,5));
t150 = sin(legFrame(5,3));
t152 = t147 * t146 + t150 * t149;
t153 = t152 * t145;
t154 = cos(legFrame(5,2));
t155 = sin(qJ(2,5));
t156 = 0.1e1 / t155;
t157 = t156 * t154;
t158 = cos(legFrame(5,1));
t160 = sin(legFrame(5,1));
t162 = -t158 * g(2) - g(3) * t160;
t164 = sin(legFrame(5,2));
t170 = t160 * g(2) * t164 - g(3) * t158 * t164 + t154 * g(1);
t172 = t170 * t146 + t149 * t162;
t175 = -t162 * t146;
t176 = t149 * t170 + t175;
t178 = t147 * t172 + t150 * t176;
t181 = 0.1e1 / qJ(3,6);
t182 = sin(qJ(1,6));
t183 = cos(legFrame(6,3));
t185 = cos(qJ(1,6));
t186 = sin(legFrame(6,3));
t188 = t183 * t182 + t186 * t185;
t189 = t188 * t181;
t190 = cos(legFrame(6,2));
t191 = sin(qJ(2,6));
t192 = 0.1e1 / t191;
t193 = t192 * t190;
t194 = cos(legFrame(6,1));
t196 = sin(legFrame(6,1));
t198 = -t194 * g(2) - g(3) * t196;
t200 = sin(legFrame(6,2));
t206 = t196 * g(2) * t200 - g(3) * t194 * t200 + t190 * g(1);
t208 = t206 * t182 + t185 * t198;
t211 = -t198 * t182;
t212 = t185 * t206 + t211;
t214 = t183 * t208 + t186 * t212;
t220 = t6 * t28;
t221 = t3 * t32 - t220;
t225 = t42 * t64;
t226 = t39 * t68 - t225;
t230 = t78 * t100;
t231 = t75 * t104 - t230;
t235 = t114 * t136;
t236 = t111 * t140 - t235;
t240 = t150 * t172;
t241 = t147 * t176 - t240;
t245 = t186 * t208;
t246 = t183 * t212 - t245;
t252 = cos(qJ(2,1));
t254 = t34 * t252 * t12;
t258 = -t2 * t6 + t5 * t3;
t262 = t10 * t258 * t252 + t20 * t11;
t263 = t1 * t262;
t266 = t2 * g(2) - t5 * t23;
t268 = t20 * t2;
t272 = (t5 * g(2) + g(3) * t268) * t6;
t277 = g(3) * t2 + t5 * t21;
t279 = t5 * t25;
t289 = t6 * (t16 * (g(2) * t268 - t5 * g(3)) + t10 * t2 * g(1));
t290 = t14 * (t3 * t266 + t272) + t3 * (t16 * t277 + t279) - t289;
t293 = t10 * t16;
t295 = t10 * t14;
t297 = g(1) * t20 - g(2) * t293 + g(3) * t295;
t298 = t252 * t297;
t299 = t11 * t290 - t298;
t302 = cos(qJ(2,2));
t304 = t70 * t302 * t48;
t308 = -t38 * t42 + t41 * t39;
t312 = t46 * t308 * t302 + t56 * t47;
t313 = t37 * t312;
t317 = -t41 * g(3) * t56 + g(2) * t38;
t319 = t56 * t38;
t323 = (t41 * g(2) + g(3) * t319) * t42;
t328 = g(3) * t38 + t41 * t57;
t330 = t41 * t61;
t340 = t42 * (t52 * (g(2) * t319 - t41 * g(3)) + t46 * t38 * g(1));
t341 = t50 * (t39 * t317 + t323) + t39 * (t52 * t328 + t330) - t340;
t344 = t46 * t52;
t346 = t46 * t50;
t348 = g(1) * t56 - g(2) * t344 + g(3) * t346;
t349 = t302 * t348;
t350 = t47 * t341 - t349;
t353 = cos(qJ(2,3));
t355 = t106 * t353 * t84;
t359 = -t74 * t78 + t77 * t75;
t363 = t82 * t359 * t353 + t92 * t83;
t364 = t73 * t363;
t367 = g(2) * t74 - t77 * t95;
t369 = t92 * t74;
t373 = (t77 * g(2) + g(3) * t369) * t78;
t378 = g(3) * t74 + t77 * t93;
t380 = t77 * t97;
t390 = t78 * (t88 * (g(2) * t369 - t77 * g(3)) + t82 * t74 * g(1));
t391 = t86 * (t75 * t367 + t373) + t75 * (t88 * t378 + t380) - t390;
t394 = t82 * t88;
t396 = t82 * t86;
t398 = g(1) * t92 - g(2) * t394 + g(3) * t396;
t399 = t353 * t398;
t400 = t83 * t391 - t399;
t403 = cos(qJ(2,4));
t405 = t142 * t403 * t120;
t409 = -t110 * t114 + t113 * t111;
t413 = t118 * t409 * t403 + t128 * t119;
t414 = t109 * t413;
t415 = t128 * t113;
t418 = g(2) * t110 - g(3) * t415;
t420 = t128 * t110;
t424 = (g(2) * t113 + g(3) * t420) * t114;
t429 = g(2) * t415 + g(3) * t110;
t432 = t118 * t113 * g(1);
t442 = t114 * (t124 * (g(2) * t420 - g(3) * t113) + t118 * t110 * g(1));
t443 = t122 * (t111 * t418 + t424) + t111 * (t124 * t429 + t432) - t442;
t446 = t118 * t124;
t448 = t118 * t122;
t450 = g(1) * t128 - g(2) * t446 + g(3) * t448;
t451 = t403 * t450;
t452 = t119 * t443 - t451;
t455 = cos(qJ(2,5));
t457 = t178 * t455 * t156;
t461 = -t146 * t150 + t149 * t147;
t465 = t154 * t461 * t455 + t164 * t155;
t466 = t145 * t465;
t467 = t164 * t149;
t470 = g(2) * t146 - g(3) * t467;
t472 = t164 * t146;
t476 = (g(2) * t149 + g(3) * t472) * t150;
t481 = g(2) * t467 + g(3) * t146;
t484 = t154 * t149 * g(1);
t494 = t150 * (t160 * (g(2) * t472 - g(3) * t149) + t154 * t146 * g(1));
t495 = t158 * (t147 * t470 + t476) + t147 * (t160 * t481 + t484) - t494;
t498 = t154 * t160;
t500 = t154 * t158;
t502 = g(1) * t164 - g(2) * t498 + g(3) * t500;
t503 = t455 * t502;
t504 = t155 * t495 - t503;
t507 = cos(qJ(2,6));
t509 = t214 * t507 * t192;
t513 = -t182 * t186 + t185 * t183;
t517 = t190 * t513 * t507 + t200 * t191;
t518 = t181 * t517;
t519 = t200 * t185;
t522 = g(2) * t182 - g(3) * t519;
t524 = t200 * t182;
t528 = (g(2) * t185 + g(3) * t524) * t186;
t533 = g(2) * t519 + g(3) * t182;
t536 = t190 * t185 * g(1);
t546 = t186 * (t196 * (g(2) * t524 - g(3) * t185) + t190 * t182 * g(1));
t547 = t194 * (t183 * t522 + t528) + t183 * (t196 * t533 + t536) - t546;
t550 = t190 * t196;
t552 = t190 * t194;
t554 = g(1) * t200 - g(2) * t550 + g(3) * t552;
t555 = t507 * t554;
t556 = t191 * t547 - t555;
t558 = -t254 * t10 * t9 - t405 * t118 * t117 - t457 * t154 * t153 - t509 * t190 * t189 - t304 * t46 * t45 - t355 * t82 * t81 + t299 * t263 + t350 * t313 + t400 * t364 + t452 * t414 + t504 * t466 + t556 * t518;
t561 = t34 * t10 * t9;
t563 = t11 * t297;
t564 = t252 * t290 + t563;
t567 = t70 * t46 * t45;
t569 = t47 * t348;
t570 = t302 * t341 + t569;
t573 = t106 * t82 * t81;
t575 = t83 * t398;
t576 = t353 * t391 + t575;
t579 = t142 * t118 * t117;
t581 = t119 * t450;
t582 = t403 * t443 + t581;
t585 = t178 * t154 * t153;
t587 = t155 * t502;
t588 = t455 * t495 + t587;
t591 = t214 * t190 * t189;
t593 = t191 * t554;
t594 = t507 * t547 + t593;
t596 = t564 * t263 + t570 * t313 + t576 * t364 + t582 * t414 + t588 * t466 + t594 * t518 + t561 + t567 + t573 + t579 + t585 + t591;
t602 = t3 * (-t5 * t26 - t31) + t220;
t608 = t39 * (-t41 * t62 - t67) + t225;
t614 = t75 * (-t77 * t98 - t103) + t230;
t620 = t111 * (-t113 * t134 - t139) + t235;
t626 = t147 * (-t149 * t170 - t175) + t240;
t632 = t183 * (-t185 * t206 - t211) + t245;
t643 = t14 * (-t3 * t266 - t272) + t3 * (-t16 * t277 - t279) + t289;
t645 = t252 * t643 - t563;
t653 = t50 * (-t39 * t317 - t323) + t39 * (-t52 * t328 - t330) + t340;
t655 = t302 * t653 - t569;
t663 = t86 * (-t75 * t367 - t373) + t75 * (-t88 * t378 - t380) + t390;
t665 = t353 * t663 - t575;
t673 = t122 * (-t111 * t418 - t424) + t111 * (-t124 * t429 - t432) + t442;
t675 = t403 * t673 - t581;
t683 = t158 * (-t147 * t470 - t476) + t147 * (-t160 * t481 - t484) + t494;
t685 = t455 * t683 - t587;
t693 = t194 * (-t183 * t522 - t528) + t183 * (-t196 * t533 - t536) + t546;
t695 = t507 * t693 - t593;
t697 = t645 * t263 + t655 * t313 + t665 * t364 + t675 * t414 + t685 * t466 + t695 * t518 - t561 - t567 - t573 - t579 - t585 - t591;
t699 = t10 * t8;
t705 = t10 * t258 * t11 - t252 * t20;
t707 = t11 * t643 + t298;
t709 = t46 * t44;
t715 = t46 * t308 * t47 - t302 * t56;
t717 = t47 * t653 + t349;
t719 = t82 * t80;
t725 = t82 * t359 * t83 - t353 * t92;
t727 = t83 * t663 + t399;
t729 = t118 * t116;
t734 = t403 * t128;
t737 = t119 * t673 + t451;
t739 = t154 * t152;
t744 = t455 * t164;
t747 = t155 * t683 + t503;
t749 = t190 * t188;
t754 = t507 * t200;
t757 = t191 * t693 + t555;
t759 = -t34 * t699 - t564 * t262 + t707 * t705 - t70 * t709 - t570 * t312 + t717 * t715 - t106 * t719 - t576 * t363 + t727 * t725 - t142 * t729 + t675 * t413 + t737 * (t118 * t409 * t119 - t734) - t178 * t739 + t685 * t465 + t747 * (t154 * t461 * t155 - t744) - t214 * t749 + t695 * t517 + t757 * (t190 * t513 * t191 - t754);
t761 = cos(xP(5));
t762 = cos(xP(6));
t763 = t762 * t761;
t765 = sin(xP(4));
t766 = sin(xP(5));
t767 = t766 * t765;
t769 = cos(xP(4));
t770 = sin(xP(6));
t772 = t762 * t767 + t770 * t769;
t774 = t766 * t769;
t777 = -t762 * t774 + t770 * t765;
t779 = -g(1) * t763 - g(2) * t772 - g(3) * t777;
t781 = t770 * t761;
t785 = t762 * t769 - t770 * t767;
t789 = t762 * t765 + t770 * t774;
t791 = g(1) * t781 - g(2) * t785 - g(3) * t789;
t794 = t765 * t761;
t796 = t769 * t761;
t798 = -g(1) * t766 + g(2) * t794 - g(3) * t796;
t804 = t16 * t20;
t806 = t14 * t258 - t8 * t804;
t807 = t1 * t806;
t808 = t34 * t12;
t811 = t52 * t56;
t813 = t50 * t308 - t44 * t811;
t814 = t37 * t813;
t815 = t70 * t48;
t818 = t88 * t92;
t820 = t86 * t359 - t80 * t818;
t821 = t73 * t820;
t822 = t106 * t84;
t827 = -t116 * t124 * t128 + t122 * t409;
t828 = t109 * t827;
t829 = t142 * t120;
t834 = -t152 * t160 * t164 + t158 * t461;
t835 = t145 * t834;
t836 = t178 * t156;
t841 = -t188 * t196 * t200 + t194 * t513;
t842 = t181 * t841;
t843 = t214 * t192;
t847 = t221 * t12;
t849 = t226 * t48;
t851 = t231 * t84;
t853 = t236 * t120;
t855 = t241 * t156;
t857 = t246 * t192;
t862 = t258 * t20;
t865 = t8 * t14 + t16 * t862;
t869 = -t10 * t11 * t16 + t252 * t865;
t870 = t1 * t869;
t873 = t308 * t56;
t876 = t44 * t50 + t52 * t873;
t880 = -t46 * t47 * t52 + t302 * t876;
t881 = t37 * t880;
t884 = t359 * t92;
t887 = t80 * t86 + t88 * t884;
t891 = -t82 * t83 * t88 + t353 * t887;
t892 = t73 * t891;
t895 = t409 * t128;
t898 = t116 * t122 + t124 * t895;
t900 = t119 * t124;
t902 = -t118 * t900 + t403 * t898;
t903 = t109 * t902;
t906 = t461 * t164;
t909 = t152 * t158 + t160 * t906;
t911 = t155 * t160;
t913 = -t154 * t911 + t455 * t909;
t914 = t145 * t913;
t917 = t513 * t200;
t920 = t188 * t194 + t196 * t917;
t922 = t191 * t196;
t924 = -t190 * t922 + t507 * t920;
t925 = t181 * t924;
t927 = t254 * t807 + t299 * t870 + t304 * t814 + t350 * t881 + t355 * t821 + t400 * t892 + t405 * t828 + t452 * t903 + t457 * t835 + t504 * t914 + t509 * t842 + t556 * t925;
t929 = t34 * t807;
t931 = t70 * t814;
t933 = t106 * t821;
t935 = t142 * t828;
t937 = t178 * t835;
t939 = t214 * t842;
t941 = t564 * t870 + t570 * t881 + t576 * t892 + t582 * t903 + t588 * t914 + t594 * t925 - t929 - t931 - t933 - t935 - t937 - t939;
t944 = t602 * t12;
t946 = t608 * t48;
t948 = t614 * t84;
t950 = t620 * t120;
t952 = t626 * t156;
t954 = t632 * t192;
t964 = t645 * t870 + t655 * t881 + t665 * t892 + t675 * t903 + t685 * t914 + t695 * t925 + t929 + t931 + t933 + t935 + t937 + t939;
t1002 = t34 * t806 - t564 * t869 + t707 * (t11 * t865 + t252 * t293) + t70 * t813 - t570 * t880 + t717 * (t302 * t344 + t47 * t876) + t106 * t820 - t576 * t891 + t727 * (t353 * t394 + t83 * t887) + t142 * t827 + t675 * t902 + t737 * (t119 * t898 + t403 * t446) + t178 * t834 + t685 * t913 + t747 * (t155 * t909 + t455 * t498) + t214 * t841 + t695 * t924 + t757 * (t191 * t920 + t507 * t550);
t1013 = t14 * t8 * t20 + t258 * t16;
t1014 = t1 * t1013;
t1019 = t50 * t44 * t56 + t308 * t52;
t1020 = t37 * t1019;
t1025 = t86 * t80 * t92 + t359 * t88;
t1026 = t73 * t1025;
t1031 = t122 * t116 * t128 + t409 * t124;
t1032 = t109 * t1031;
t1037 = t158 * t152 * t164 + t461 * t160;
t1038 = t145 * t1037;
t1043 = t194 * t188 * t200 + t513 * t196;
t1044 = t181 * t1043;
t1059 = -t14 * t862 + t8 * t16;
t1061 = t11 * t14;
t1063 = t10 * t1061 + t252 * t1059;
t1064 = t1 * t1063;
t1069 = t44 * t52 - t50 * t873;
t1071 = t47 * t50;
t1073 = t302 * t1069 + t46 * t1071;
t1074 = t37 * t1073;
t1079 = t80 * t88 - t86 * t884;
t1081 = t83 * t86;
t1083 = t353 * t1079 + t82 * t1081;
t1084 = t73 * t1083;
t1089 = t116 * t124 - t122 * t895;
t1091 = t119 * t122;
t1093 = t403 * t1089 + t118 * t1091;
t1094 = t109 * t1093;
t1099 = t152 * t160 - t158 * t906;
t1101 = t155 * t158;
t1103 = t455 * t1099 + t154 * t1101;
t1104 = t145 * t1103;
t1109 = t188 * t196 - t194 * t917;
t1111 = t191 * t194;
t1113 = t507 * t1109 + t190 * t1111;
t1114 = t181 * t1113;
t1116 = t254 * t1014 + t304 * t1020 + t355 * t1026 + t405 * t1032 + t457 * t1038 + t509 * t1044 + t299 * t1064 + t350 * t1074 + t400 * t1084 + t452 * t1094 + t504 * t1104 + t556 * t1114;
t1118 = t34 * t1014;
t1120 = t70 * t1020;
t1122 = t106 * t1026;
t1124 = t142 * t1032;
t1126 = t178 * t1038;
t1128 = t214 * t1044;
t1130 = t564 * t1064 + t570 * t1074 + t576 * t1084 + t582 * t1094 + t588 * t1104 + t594 * t1114 - t1118 - t1120 - t1122 - t1124 - t1126 - t1128;
t1147 = t645 * t1064 + t655 * t1074 + t665 * t1084 + t675 * t1094 + t685 * t1104 + t695 * t1114 + t1118 + t1120 + t1122 + t1124 + t1126 + t1128;
t1185 = t34 * t1013 - t564 * t1063 + t707 * (t11 * t1059 - t252 * t295) + t70 * t1019 - t570 * t1073 + t717 * (t47 * t1069 - t302 * t346) + t106 * t1025 - t576 * t1083 + t727 * (t83 * t1079 - t353 * t396) + t142 * t1031 + t675 * t1093 + t737 * (t119 * t1089 - t403 * t448) + t178 * t1037 + t685 * t1103 + t747 * (t155 * t1099 - t455 * t500) + t214 * t1043 + t695 * t1113 + t757 * (t191 * t1109 - t507 * t552);
t1193 = t769 * koppelP(1,2);
t1195 = t765 * koppelP(1,1);
t1196 = -t766 * t1193 - t1195;
t1198 = t769 * koppelP(1,1);
t1200 = t765 * koppelP(1,2);
t1201 = t766 * t1198 - t1200;
t1204 = t761 * t769 * koppelP(1,3);
t1205 = t770 * t1196 + t762 * t1201 - t1204;
t1208 = -t14 * t5 + t16 * t268;
t1212 = t14 * t2 + t5 * t804;
t1216 = t12 * t1;
t1218 = t14 * t20;
t1221 = t2 * t1218 + t16 * t5;
t1225 = t5 * t1218 - t16 * t2;
t1234 = t761 * koppelP(1,3);
t1236 = t770 * (-t766 * t1200 + t1198) + t762 * (t766 * t1195 + t1193) - t765 * t1234;
t1239 = -t1216 * (t3 * t1208 + t1212 * t6) * t1205 + t1216 * t1236 * (t3 * t1221 + t1225 * t6);
t1242 = koppelP(2,2) * t765;
t1243 = koppelP(2,1) * t774 - t1242;
t1247 = -koppelP(2,1) * t765 - koppelP(2,2) * t774;
t1249 = koppelP(2,3) * t796;
t1250 = t762 * t1243 + t770 * t1247 - t1249;
t1253 = t52 * t319 - t41 * t50;
t1257 = t50 * t38 + t41 * t811;
t1261 = t48 * t37;
t1265 = t50 * t319 + t41 * t52;
t1268 = t52 * t38;
t1269 = t41 * t59 - t1268;
t1282 = -t761 * koppelP(2,3) * t765 + t762 * (koppelP(2,1) * t767 + koppelP(2,2) * t769) + (-t766 * t1242 + koppelP(2,1) * t769) * t770;
t1285 = -t1261 * (t39 * t1253 + t1257 * t42) * t1250 + t1261 * t1282 * (t39 * t1265 + t1269 * t42);
t1288 = koppelP(3,2) * t765;
t1289 = koppelP(3,1) * t774 - t1288;
t1293 = -koppelP(3,1) * t765 - koppelP(3,2) * t774;
t1295 = koppelP(3,3) * t796;
t1296 = t762 * t1289 + t770 * t1293 - t1295;
t1299 = t88 * t369 - t77 * t86;
t1303 = t86 * t74 + t77 * t818;
t1307 = t84 * t73;
t1311 = t86 * t369 + t77 * t88;
t1315 = t74 * t88;
t1316 = t77 * t86 * t92 - t1315;
t1329 = -t761 * koppelP(3,3) * t765 + t762 * (koppelP(3,1) * t767 + koppelP(3,2) * t769) + (-t766 * t1288 + koppelP(3,1) * t769) * t770;
t1332 = -t1307 * (t75 * t1299 + t1303 * t78) * t1296 + t1307 * t1329 * (t75 * t1311 + t1316 * t78);
t1333 = t106 * t1332;
t1336 = -t122 * t113 + t124 * t420;
t1340 = t122 * t110 + t124 * t415;
t1344 = koppelP(4,2) * t765;
t1345 = koppelP(4,1) * t774 - t1344;
t1349 = -koppelP(4,1) * t765 - koppelP(4,2) * t774;
t1351 = koppelP(4,3) * t796;
t1352 = t762 * t1345 + t770 * t1349 - t1351;
t1354 = t120 * t109;
t1358 = t124 * t113 + t122 * t420;
t1362 = -t124 * t110 + t122 * t415;
t1375 = -t761 * koppelP(4,3) * t765 + t762 * (koppelP(4,1) * t767 + koppelP(4,2) * t769) + (-t766 * t1344 + koppelP(4,1) * t769) * t770;
t1378 = -t1354 * t1352 * (t111 * t1336 + t1340 * t114) + t1354 * t1375 * (t111 * t1358 + t1362 * t114);
t1382 = -t158 * t149 + t160 * t472;
t1386 = t158 * t146 + t160 * t467;
t1390 = koppelP(5,2) * t765;
t1391 = koppelP(5,1) * t774 - t1390;
t1395 = -koppelP(5,1) * t765 - koppelP(5,2) * t774;
t1397 = koppelP(5,3) * t796;
t1398 = t762 * t1391 + t770 * t1395 - t1397;
t1400 = t156 * t145;
t1412 = -t761 * koppelP(5,3) * t765 + t762 * (koppelP(5,1) * t767 + koppelP(5,2) * t769) + (-t766 * t1390 + koppelP(5,1) * t769) * t770;
t1415 = t160 * t149 + t158 * t472;
t1419 = -t160 * t146 + t158 * t467;
t1424 = -t1400 * t1398 * (t147 * t1382 + t1386 * t150) + t1400 * (t147 * t1415 + t1419 * t150) * t1412;
t1427 = koppelP(6,2) * t765;
t1428 = koppelP(6,1) * t774 - t1427;
t1432 = -koppelP(6,1) * t765 - koppelP(6,2) * t774;
t1434 = koppelP(6,3) * t796;
t1435 = t762 * t1428 + t770 * t1432 - t1434;
t1438 = -t194 * t185 + t196 * t524;
t1442 = t194 * t182 + t196 * t519;
t1446 = t192 * t181;
t1458 = -t761 * koppelP(6,3) * t765 + t762 * (koppelP(6,1) * t767 + koppelP(6,2) * t769) + (-t766 * t1427 + koppelP(6,1) * t769) * t770;
t1461 = t196 * t185 + t194 * t524;
t1465 = -t196 * t182 + t194 * t519;
t1470 = -t1446 * (t183 * t1438 + t1442 * t186) * t1435 + t1446 * (t183 * t1461 + t1465 * t186) * t1458;
t1486 = t762 * koppelP(1,2) + koppelP(1,1) * t770;
t1487 = t1486 * t5;
t1491 = t770 * t766;
t1493 = -t762 * koppelP(1,1) * t766 + koppelP(1,2) * t1491 + t1234;
t1494 = t1493 * t2;
t1495 = -t20 * t1487 - t1494;
t1497 = t1486 * t2;
t1499 = t1493 * t5;
t1501 = t6 * (-t20 * t1497 + t1499);
t1510 = t3 * (t20 * t1499 - t1497) - (t20 * t1494 + t1487) * t6;
t1533 = t252 * (t14 * (t769 * (t3 * t1495 - t1501) + t1510 * t765) - t16 * (t769 * t1510 + (-t3 * t1495 + t1501) * t765)) + (t14 * (t1486 * t769 - t1493 * t765) + (t1486 * t765 + t769 * t1493) * t16) * t11 * t10;
t1534 = t1 * t1533;
t1540 = t762 * koppelP(2,2) + koppelP(2,1) * t770;
t1541 = t1540 * t41;
t1547 = -t762 * koppelP(2,1) * t766 + koppelP(2,2) * t1491 + koppelP(2,3) * t761;
t1548 = t1547 * t38;
t1549 = -t56 * t1541 - t1548;
t1551 = t1540 * t38;
t1553 = t1547 * t41;
t1555 = (-t56 * t1551 + t1553) * t42;
t1564 = t39 * (t56 * t1553 - t1551) - (t56 * t1548 + t1541) * t42;
t1587 = t302 * (t50 * (t769 * (t39 * t1549 - t1555) + t1564 * t765) - t52 * (t769 * t1564 + (-t39 * t1549 + t1555) * t765)) + (t50 * (t1540 * t769 - t1547 * t765) + (t1540 * t765 + t769 * t1547) * t52) * t47 * t46;
t1588 = t37 * t1587;
t1594 = t762 * koppelP(3,2) + koppelP(3,1) * t770;
t1595 = t1594 * t77;
t1601 = -t762 * koppelP(3,1) * t766 + koppelP(3,2) * t1491 + koppelP(3,3) * t761;
t1602 = t1601 * t74;
t1603 = -t92 * t1595 - t1602;
t1605 = t1594 * t74;
t1607 = t1601 * t77;
t1609 = (-t92 * t1605 + t1607) * t78;
t1618 = t75 * (t92 * t1607 - t1605) - (t92 * t1602 + t1595) * t78;
t1641 = t353 * (t86 * (t769 * (t75 * t1603 - t1609) + t1618 * t765) - t88 * (t769 * t1618 + (-t75 * t1603 + t1609) * t765)) + (t86 * (t1594 * t769 - t1601 * t765) + (t1594 * t765 + t769 * t1601) * t88) * t83 * t82;
t1642 = t73 * t1641;
t1648 = t762 * koppelP(4,2) + koppelP(4,1) * t770;
t1649 = t1648 * t113;
t1655 = -t762 * koppelP(4,1) * t766 + koppelP(4,2) * t1491 + koppelP(4,3) * t761;
t1656 = t1655 * t110;
t1657 = -t128 * t1649 - t1656;
t1659 = t1648 * t110;
t1661 = t1655 * t113;
t1663 = (-t128 * t1659 + t1661) * t114;
t1672 = t111 * (t128 * t1661 - t1659) - (t128 * t1656 + t1649) * t114;
t1695 = t403 * (t122 * (t769 * (t111 * t1657 - t1663) + t765 * t1672) - t124 * (t769 * t1672 + t765 * (-t111 * t1657 + t1663))) + (t122 * (t1648 * t769 - t1655 * t765) + (t1648 * t765 + t769 * t1655) * t124) * t119 * t118;
t1696 = t109 * t1695;
t1702 = t762 * koppelP(5,2) + koppelP(5,1) * t770;
t1703 = t1702 * t149;
t1709 = -t762 * koppelP(5,1) * t766 + koppelP(5,2) * t1491 + koppelP(5,3) * t761;
t1710 = t1709 * t146;
t1711 = -t164 * t1703 - t1710;
t1713 = t1702 * t146;
t1715 = t1709 * t149;
t1717 = (-t164 * t1713 + t1715) * t150;
t1726 = t147 * (t164 * t1715 - t1713) - (t164 * t1710 + t1703) * t150;
t1749 = t455 * (t158 * (t769 * (t147 * t1711 - t1717) + t765 * t1726) - t160 * (t769 * t1726 + t765 * (-t147 * t1711 + t1717))) + (t158 * (t1702 * t769 - t1709 * t765) + (t1702 * t765 + t769 * t1709) * t160) * t155 * t154;
t1750 = t145 * t1749;
t1756 = t762 * koppelP(6,2) + koppelP(6,1) * t770;
t1757 = t1756 * t185;
t1763 = -t762 * koppelP(6,1) * t766 + koppelP(6,2) * t1491 + koppelP(6,3) * t761;
t1764 = t1763 * t182;
t1765 = -t200 * t1757 - t1764;
t1767 = t1756 * t182;
t1769 = t1763 * t185;
t1771 = (-t200 * t1767 + t1769) * t186;
t1780 = t183 * (t200 * t1769 - t1767) - (t200 * t1764 + t1757) * t186;
t1803 = t507 * (t194 * (t769 * (t183 * t1765 - t1771) + t765 * t1780) - t196 * (t769 * t1780 + t765 * (-t183 * t1765 + t1771))) + (t194 * (t1756 * t769 - t1763 * t765) + (t1756 * t765 + t769 * t1763) * t196) * t191 * t190;
t1804 = t181 * t1803;
t1806 = t106 * t353 * t1332 + t34 * t252 * t1239 + t70 * t302 * t1285 + t142 * t403 * t1378 + t178 * t455 * t1424 + t214 * t507 * t1470 + t299 * t1534 + t350 * t1588 + t400 * t1642 + t452 * t1696 + t504 * t1750 + t556 * t1804;
t1808 = t11 * t1239;
t1809 = t34 * t1808;
t1811 = t47 * t1285;
t1812 = t70 * t1811;
t1814 = t83 * t1333;
t1817 = t142 * t119 * t1378;
t1820 = t178 * t155 * t1424;
t1823 = t214 * t191 * t1470;
t1825 = t564 * t1534 + t570 * t1588 + t576 * t1642 + t582 * t1696 + t588 * t1750 + t594 * t1804 - t1809 - t1812 - t1814 - t1817 - t1820 - t1823;
t1842 = t645 * t1534 + t655 * t1588 + t665 * t1642 + t675 * t1696 + t685 * t1750 + t695 * t1804 + t1809 + t1812 + t1814 + t1817 + t1820 + t1823;
t1844 = qJ(3,1) * t34;
t1849 = t3 * t20 * t11 * t5;
t1850 = t11 * t2;
t1853 = t10 * t252;
t1857 = t16 * (-t6 * t20 * t1850 + t1849 + t1853) + t8 * t1061;
t1865 = -t6 * t1221 * t11 - t1850 * t16 * t3 + (t1849 + t1853) * t14;
t1869 = qJ(3,2) * t70;
t1874 = t39 * t56 * t47 * t41;
t1878 = t46 * t302;
t1882 = t52 * (-t42 * t56 * t47 * t38 + t1874 + t1878) + t44 * t1071;
t1890 = -t42 * t1265 * t47 - t1268 * t39 * t47 + (t1874 + t1878) * t50;
t1894 = qJ(3,3) * t83;
t1899 = t75 * t92 * t83 * t77;
t1903 = t82 * t353;
t1907 = t88 * (-t78 * t92 * t83 * t74 + t1899 + t1903) + t80 * t1081;
t1915 = -t78 * t1311 * t83 - t1315 * t75 * t83 + (t1899 + t1903) * t86;
t1920 = t142 * t119;
t1925 = t119 * t110;
t1927 = t114 * t128 * t1925;
t1928 = t118 * t403;
t1931 = t119 * t113;
t1938 = t111 * t128 * t1931 - t1927 + t1928;
t1941 = -t116 * t900 + t122 * t1938;
t1946 = t178 * t155;
t1951 = t155 * t146;
t1953 = t150 * t164 * t1951;
t1954 = t154 * t455;
t1957 = t155 * t149;
t1964 = t147 * t164 * t1957 - t1953 + t1954;
t1967 = -t152 * t911 + t158 * t1964;
t1972 = t214 * t191;
t1979 = t186 * t200 * t191 * t182;
t1980 = t190 * t507;
t1983 = t191 * t185;
t1990 = t183 * t200 * t1983 - t1979 + t1980;
t1993 = -t188 * t922 + t194 * t1990;
t1997 = t1844 * t1808 - t564 * t1533 + t707 * (t1857 * t1205 - t1236 * t1865) + t1869 * t1811 - t570 * t1587 + t717 * (t1882 * t1250 - t1282 * t1890) + t1894 * t1333 - t576 * t1641 + t727 * (t1907 * t1296 - t1329 * t1915) + t1920 * qJ(3,4) * t1378 + t675 * t1695 + t737 * ((t111 * t1340 * t119 + t124 * (-t1927 + t1928) + t114 * t122 * t1931) * t1352 - t1375 * t1941) + t1946 * qJ(3,5) * t1424 + t685 * t1749 + t747 * ((t147 * t1386 * t155 + t160 * (-t1953 + t1954) + t150 * t158 * t1957) * t1398 - t1412 * t1967) + t1972 * qJ(3,6) * t1470 + t695 * t1803 + t757 * ((t183 * t1442 * t191 + t196 * (-t1979 + t1980) + t186 * t194 * t1983) * t1435 - t1458 * t1993);
t2011 = MDP(2) * (t34 * t1239 + t70 * t1285 + t142 * t1378 + t178 * t1424 + t214 * t1470 + t1333) + MDP(3) * (t221 * t1239 + t226 * t1285 + t231 * t1332 + t236 * t1378 + t241 * t1424 + t246 * t1470) + MDP(9) * t1806 + MDP(10) * t1825 + MDP(11) * t1806 + MDP(12) * (t602 * t1239 + t608 * t1285 + t614 * t1332 + t620 * t1378 + t626 * t1424 + t632 * t1470) + MDP(13) * t1842 + MDP(14) * t1997 + MDP(21) * (t791 * t766 + t798 * t781) + MDP(22) * (t798 * t763 - t779 * t766) + MDP(23) * (-t791 * t763 - t779 * t781);
t2017 = koppelP(1,3) * t766;
t2018 = -koppelP(1,1) * t763 + koppelP(1,2) * t781 - t2017;
t2021 = t12 * t1 * t1205 * t699 + t1216 * t1013 * t2018;
t2028 = koppelP(2,3) * t766;
t2029 = -koppelP(2,1) * t763 + koppelP(2,2) * t781 - t2028;
t2032 = t48 * t37 * t1250 * t709 + t1261 * t1019 * t2029;
t2039 = koppelP(3,3) * t766;
t2040 = -koppelP(3,1) * t763 + koppelP(3,2) * t781 - t2039;
t2043 = t84 * t73 * t1296 * t719 + t1307 * t1025 * t2040;
t2044 = t106 * t2043;
t2050 = koppelP(4,3) * t766;
t2051 = -koppelP(4,1) * t763 + koppelP(4,2) * t781 - t2050;
t2054 = t120 * t109 * t1352 * t729 + t1354 * t1031 * t2051;
t2061 = koppelP(5,3) * t766;
t2062 = -koppelP(5,1) * t763 + koppelP(5,2) * t781 - t2061;
t2065 = t156 * t145 * t1398 * t739 + t1400 * t1037 * t2062;
t2072 = koppelP(6,3) * t766;
t2073 = -koppelP(6,1) * t763 + koppelP(6,2) * t781 - t2072;
t2076 = t192 * t181 * t1435 * t749 + t1446 * t1043 * t2073;
t2092 = -t770 * t1196 - t762 * t1201 + t1204;
t2099 = t761 * (t762 * koppelP(1,1) - koppelP(1,2) * t770) + t2017;
t2111 = t252 * (t10 * t258 * t2092 + (-t1221 * t6 + t3 * t1225) * t2099) + t11 * (-t10 * t14 * t2099 + t2092 * t20);
t2112 = t1 * t2111;
t2118 = -t762 * t1243 - t770 * t1247 + t1249;
t2125 = t761 * (t762 * koppelP(2,1) - koppelP(2,2) * t770) + t2028;
t2137 = t302 * (t46 * t308 * t2118 + (-t1265 * t42 + t39 * t1269) * t2125) + t47 * (-t46 * t50 * t2125 + t2118 * t56);
t2138 = t37 * t2137;
t2144 = -t762 * t1289 - t770 * t1293 + t1295;
t2151 = t761 * (t762 * koppelP(3,1) - koppelP(3,2) * t770) + t2039;
t2163 = t353 * (t82 * t359 * t2144 + (-t1311 * t78 + t75 * t1316) * t2151) + t83 * (-t82 * t86 * t2151 + t2144 * t92);
t2164 = t73 * t2163;
t2170 = -t762 * t1345 - t770 * t1349 + t1351;
t2177 = t761 * (t762 * koppelP(4,1) - koppelP(4,2) * t770) + t2050;
t2189 = t403 * (t118 * t2170 * t409 + (t111 * t1362 - t1358 * t114) * t2177) + (-t118 * t122 * t2177 + t2170 * t128) * t119;
t2190 = t109 * t2189;
t2196 = -t762 * t1391 - t770 * t1395 + t1397;
t2203 = t761 * (t762 * koppelP(5,1) - koppelP(5,2) * t770) + t2061;
t2215 = t455 * (t154 * t2196 * t461 + (-t1415 * t150 + t147 * t1419) * t2203) + (-t154 * t158 * t2203 + t2196 * t164) * t155;
t2216 = t145 * t2215;
t2222 = -t762 * t1428 - t770 * t1432 + t1434;
t2229 = t761 * (t762 * koppelP(6,1) - koppelP(6,2) * t770) + t2072;
t2241 = t507 * (t190 * t2222 * t513 + (-t1461 * t186 + t183 * t1465) * t2229) + (-t190 * t194 * t2229 + t2222 * t200) * t191;
t2242 = t181 * t2241;
t2244 = t106 * t353 * t2043 + t142 * t403 * t2054 + t178 * t455 * t2065 + t34 * t252 * t2021 + t70 * t302 * t2032 + t214 * t507 * t2076 + t299 * t2112 + t350 * t2138 + t400 * t2164 + t452 * t2190 + t504 * t2216 + t556 * t2242;
t2246 = t11 * t2021;
t2247 = t34 * t2246;
t2249 = t47 * t2032;
t2250 = t70 * t2249;
t2252 = t83 * t2044;
t2255 = t142 * t119 * t2054;
t2258 = t178 * t155 * t2065;
t2261 = t214 * t191 * t2076;
t2263 = t564 * t2112 + t570 * t2138 + t576 * t2164 + t582 * t2190 + t588 * t2216 + t594 * t2242 - t2247 - t2250 - t2252 - t2255 - t2258 - t2261;
t2280 = t645 * t2112 + t655 * t2138 + t665 * t2164 + t675 * t2190 + t685 * t2216 + t695 * t2242 + t2247 + t2250 + t2252 + t2255 + t2258 + t2261;
t2307 = t111 * t118 * t1931 - t114 * t118 * t1925 - t734;
t2319 = t147 * t154 * t1957 - t150 * t154 * t1951 - t744;
t2333 = -t190 * t182 * t186 * t191 + t190 * t185 * t183 * t191 - t754;
t2338 = t1844 * t2246 - t564 * t2111 + t707 * (-t705 * t1205 + t1865 * t2099) + t1869 * t2249 - t570 * t2137 + t717 * (-t715 * t1250 + t1890 * t2125) + t1894 * t2044 - t576 * t2163 + t727 * (-t725 * t1296 + t1915 * t2151) + t1920 * qJ(3,4) * t2054 + t675 * t2189 + t737 * (-t1352 * t2307 + t1941 * t2177) + t1946 * qJ(3,5) * t2065 + t685 * t2215 + t747 * (-t1398 * t2319 + t1967 * t2203) + t1972 * qJ(3,6) * t2076 + t695 * t2241 + t757 * (-t1435 * t2333 + t1993 * t2229);
t2352 = MDP(2) * (t142 * t2054 + t178 * t2065 + t34 * t2021 + t70 * t2032 + t214 * t2076 + t2044) + MDP(3) * (t221 * t2021 + t226 * t2032 + t231 * t2043 + t236 * t2054 + t241 * t2065 + t246 * t2076) + MDP(9) * t2244 + MDP(10) * t2263 + MDP(11) * t2244 + MDP(12) * (t602 * t2021 + t608 * t2032 + t614 * t2043 + t620 * t2054 + t626 * t2065 + t632 * t2076) + MDP(13) * t2280 + MDP(14) * t2338 + MDP(21) * (-t798 * t785 - t791 * t794) + MDP(22) * (t798 * t772 + t779 * t794) + MDP(23) * (-t791 * t772 + t779 * t785);
t2358 = t12 * t1 * t1236 * t699 - t1216 * t806 * t2018;
t2365 = t48 * t37 * t1282 * t709 - t1261 * t813 * t2029;
t2372 = t84 * t73 * t1329 * t719 - t1307 * t820 * t2040;
t2373 = t106 * t2372;
t2379 = t120 * t109 * t1375 * t729 - t1354 * t827 * t2051;
t2386 = t156 * t145 * t1412 * t739 - t1400 * t834 * t2062;
t2393 = t192 * t181 * t1458 * t749 - t1446 * t841 * t2073;
t2420 = t252 * (-t10 * t258 * t1236 + t2099 * (-t1208 * t6 + t3 * t1212)) - (t10 * t16 * t2099 + t1236 * t20) * t11;
t2421 = t1 * t2420;
t2438 = t302 * (-t46 * t308 * t1282 + (-t1253 * t42 + t39 * t1257) * t2125) - (t46 * t52 * t2125 + t1282 * t56) * t47;
t2439 = t37 * t2438;
t2456 = t353 * (-t82 * t359 * t1329 + (-t1299 * t78 + t75 * t1303) * t2151) - (t82 * t88 * t2151 + t1329 * t92) * t83;
t2457 = t73 * t2456;
t2474 = t403 * (-t118 * t409 * t1375 + (t111 * t1340 - t1336 * t114) * t2177) - t119 * (t118 * t124 * t2177 + t1375 * t128);
t2475 = t109 * t2474;
t2492 = t455 * (-t154 * t461 * t1412 + (-t1382 * t150 + t147 * t1386) * t2203) - t155 * (t154 * t160 * t2203 + t1412 * t164);
t2493 = t145 * t2492;
t2510 = t507 * (-t190 * t513 * t1458 + (-t1438 * t186 + t183 * t1442) * t2229) - t191 * (t190 * t196 * t2229 + t1458 * t200);
t2511 = t181 * t2510;
t2513 = t106 * t353 * t2372 + t142 * t403 * t2379 + t178 * t455 * t2386 + t214 * t507 * t2393 + t34 * t252 * t2358 + t70 * t302 * t2365 + t299 * t2421 + t350 * t2439 + t400 * t2457 + t452 * t2475 + t504 * t2493 + t556 * t2511;
t2515 = t11 * t2358;
t2516 = t34 * t2515;
t2518 = t47 * t2365;
t2519 = t70 * t2518;
t2521 = t83 * t2373;
t2524 = t142 * t119 * t2379;
t2527 = t178 * t155 * t2386;
t2530 = t214 * t191 * t2393;
t2532 = t564 * t2421 + t570 * t2439 + t576 * t2457 + t582 * t2475 + t588 * t2493 + t594 * t2511 - t2516 - t2519 - t2521 - t2524 - t2527 - t2530;
t2549 = t645 * t2421 + t655 * t2439 + t665 * t2457 + t675 * t2475 + t685 * t2493 + t695 * t2511 + t2516 + t2519 + t2521 + t2524 + t2527 + t2530;
t2599 = t1844 * t2515 - t564 * t2420 + t707 * (-t705 * t1236 + t1857 * t2099) + t1869 * t2518 - t570 * t2438 + t717 * (-t715 * t1282 + t1882 * t2125) + t1894 * t2373 - t576 * t2456 + t727 * (-t725 * t1329 + t1907 * t2151) + t1920 * qJ(3,4) * t2379 + t675 * t2474 + t737 * (-t2307 * t1375 + t2177 * (t116 * t1091 + t124 * t1938)) + t1946 * qJ(3,5) * t2386 + t685 * t2492 + t747 * (-t2319 * t1412 + t2203 * (t152 * t1101 + t160 * t1964)) + t1972 * qJ(3,6) * t2393 + t695 * t2510 + t757 * (-t2333 * t1458 + t2229 * (t188 * t1111 + t196 * t1990));
t2613 = MDP(2) * (t142 * t2379 + t178 * t2386 + t214 * t2393 + t34 * t2358 + t70 * t2365 + t2373) + MDP(3) * (t221 * t2358 + t226 * t2365 + t231 * t2372 + t236 * t2379 + t241 * t2386 + t246 * t2393) + MDP(9) * t2513 + MDP(10) * t2532 + MDP(11) * t2513 + MDP(12) * (t602 * t2358 + t608 * t2365 + t614 * t2372 + t620 * t2379 + t626 * t2386 + t632 * t2393) + MDP(13) * t2549 + MDP(14) * t2599 + MDP(21) * (-t798 * t789 + t791 * t796) + MDP(22) * (t798 * t777 - t779 * t796) + MDP(23) * (-t791 * t777 + t779 * t789);
unknown(1,1) = MDP(2) * (-t106 * t85 * t81 - t142 * t121 * t117 - t34 * t13 * t9 - t178 * t157 * t153 - t214 * t193 * t189 - t70 * t49 * t45) + MDP(3) * (-t236 * t121 * t117 - t221 * t13 * t9 - t241 * t157 * t153 - t246 * t193 * t189 - t226 * t49 * t45 - t231 * t85 * t81) + MDP(9) * t558 + MDP(10) * t596 + MDP(11) * t558 + MDP(12) * (-t620 * t121 * t117 - t602 * t13 * t9 - t626 * t157 * t153 - t632 * t193 * t189 - t608 * t49 * t45 - t614 * t85 * t81) + MDP(13) * t697 + MDP(14) * t759 + MDP(24) * (t779 * t763 + t798 * t766 - t791 * t781);
unknown(2,1) = MDP(2) * (t808 * t807 + t815 * t814 + t822 * t821 + t829 * t828 + t836 * t835 + t843 * t842) + MDP(3) * (t847 * t807 + t849 * t814 + t851 * t821 + t853 * t828 + t855 * t835 + t857 * t842) + MDP(9) * t927 + MDP(10) * t941 + MDP(11) * t927 + MDP(12) * (t944 * t807 + t946 * t814 + t948 * t821 + t950 * t828 + t952 * t835 + t954 * t842) + MDP(13) * t964 + MDP(14) * t1002 + MDP(24) * (t779 * t772 + t791 * t785 - t798 * t794);
unknown(3,1) = MDP(2) * (t808 * t1014 + t815 * t1020 + t822 * t1026 + t829 * t1032 + t836 * t1038 + t843 * t1044) + MDP(3) * (t847 * t1014 + t849 * t1020 + t851 * t1026 + t853 * t1032 + t855 * t1038 + t857 * t1044) + MDP(9) * t1116 + MDP(10) * t1130 + MDP(11) * t1116 + MDP(12) * (t944 * t1014 + t946 * t1020 + t948 * t1026 + t950 * t1032 + t952 * t1038 + t954 * t1044) + MDP(13) * t1147 + MDP(14) * t1185 + MDP(24) * (t779 * t777 + t791 * t789 + t798 * t796);
unknown(4,1) = t2011;
unknown(5,1) = t2352;
unknown(6,1) = t2613;
taugX  = unknown;
