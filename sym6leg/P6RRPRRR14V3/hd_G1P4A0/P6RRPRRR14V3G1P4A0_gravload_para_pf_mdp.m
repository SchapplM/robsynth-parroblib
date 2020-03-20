% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P6RRPRRR14V3G1P4A0
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
%   see P6RRPRRR14V3G1P4A0_convert_par2_MPV_fixb.m

% Output:
% taugX [6x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 23:28
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(24,1)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: qJ has to be [3x6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: Koppelpunkt has to be [6x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_mdp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 1
% StartTime: 2020-03-12 23:28:44
% EndTime: 2020-03-12 23:28:44
% DurationCPUTime: 0.63s
% Computational Cost: add. (4251->679), mult. (9051->1134), div. (660->12), fcn. (9486->42), ass. (0->442)
unknown=NaN(6,1);
t1 = sin(legFrame(1,3));
t2 = cos(qJ(1,1));
t4 = sin(qJ(1,1));
t5 = cos(legFrame(1,3));
t7 = -t2 * t1 - t5 * t4;
t8 = sin(qJ(2,1));
t9 = 0.1e1 / t8;
t10 = t9 * t7;
t11 = 0.1e1 / qJ(3,1);
t14 = g(1) * t1 - g(2) * t5;
t18 = g(1) * t5 + g(2) * t1;
t20 = t2 * t14 + t4 * t18;
t21 = t20 * t11;
t23 = sin(legFrame(2,3));
t24 = cos(qJ(1,2));
t26 = sin(qJ(1,2));
t27 = cos(legFrame(2,3));
t29 = -t24 * t23 - t27 * t26;
t30 = sin(qJ(2,2));
t31 = 0.1e1 / t30;
t32 = t31 * t29;
t33 = 0.1e1 / qJ(3,2);
t36 = g(1) * t23 - g(2) * t27;
t40 = g(1) * t27 + g(2) * t23;
t42 = t24 * t36 + t26 * t40;
t43 = t42 * t33;
t45 = sin(legFrame(3,3));
t46 = cos(qJ(1,3));
t48 = sin(qJ(1,3));
t49 = cos(legFrame(3,3));
t51 = -t46 * t45 - t49 * t48;
t52 = sin(qJ(2,3));
t53 = 0.1e1 / t52;
t54 = t53 * t51;
t55 = 0.1e1 / qJ(3,3);
t58 = g(1) * t45 - g(2) * t49;
t62 = g(1) * t49 + g(2) * t45;
t64 = t46 * t58 + t48 * t62;
t65 = t64 * t55;
t67 = sin(legFrame(4,3));
t68 = cos(qJ(1,4));
t70 = sin(qJ(1,4));
t71 = cos(legFrame(4,3));
t73 = -t68 * t67 - t71 * t70;
t74 = sin(qJ(2,4));
t75 = 0.1e1 / t74;
t76 = t75 * t73;
t77 = 0.1e1 / qJ(3,4);
t80 = g(1) * t67 - g(2) * t71;
t84 = g(1) * t71 + g(2) * t67;
t86 = t68 * t80 + t70 * t84;
t87 = t86 * t77;
t89 = sin(legFrame(5,3));
t90 = cos(qJ(1,5));
t92 = sin(qJ(1,5));
t93 = cos(legFrame(5,3));
t95 = -t90 * t89 - t93 * t92;
t96 = sin(qJ(2,5));
t97 = 0.1e1 / t96;
t98 = t97 * t95;
t99 = 0.1e1 / qJ(3,5);
t102 = g(1) * t89 - g(2) * t93;
t106 = g(1) * t93 + g(2) * t89;
t108 = t90 * t102 + t92 * t106;
t109 = t108 * t99;
t111 = sin(legFrame(6,3));
t112 = cos(qJ(1,6));
t114 = sin(qJ(1,6));
t115 = cos(legFrame(6,3));
t117 = -t112 * t111 - t115 * t114;
t118 = sin(qJ(2,6));
t119 = 0.1e1 / t118;
t120 = t119 * t117;
t121 = 0.1e1 / qJ(3,6);
t124 = g(1) * t111 - g(2) * t115;
t128 = g(1) * t115 + g(2) * t111;
t130 = t112 * t124 + t114 * t128;
t131 = t130 * t121;
t135 = t2 * t18;
t137 = -t4 * t14 + t135;
t138 = t137 * t11;
t140 = t24 * t40;
t142 = -t26 * t36 + t140;
t143 = t142 * t33;
t145 = t46 * t62;
t147 = -t48 * t58 + t145;
t148 = t147 * t55;
t150 = t68 * t84;
t152 = -t70 * t80 + t150;
t153 = t152 * t77;
t155 = t90 * t106;
t157 = -t92 * t102 + t155;
t158 = t157 * t99;
t160 = t112 * t128;
t162 = -t114 * t124 + t160;
t163 = t162 * t121;
t167 = cos(qJ(2,1));
t169 = t20 * t167 * t11;
t173 = -t4 * t1 + t5 * t2;
t174 = t167 * t173;
t177 = t2 * g(1) + t4 * g(2);
t182 = (t4 * g(1) - t2 * g(2)) * t1;
t183 = t5 * t177 - t182;
t185 = t167 * g(3);
t186 = t8 * t183 - t185;
t189 = cos(qJ(2,2));
t191 = t42 * t189 * t33;
t195 = -t26 * t23 + t27 * t24;
t196 = t189 * t195;
t199 = t24 * g(1) + g(2) * t26;
t204 = (t26 * g(1) - t24 * g(2)) * t23;
t205 = t27 * t199 - t204;
t207 = t189 * g(3);
t208 = t30 * t205 - t207;
t211 = cos(qJ(2,3));
t213 = t64 * t211 * t55;
t217 = -t48 * t45 + t49 * t46;
t218 = t211 * t217;
t221 = t46 * g(1) + g(2) * t48;
t226 = (t48 * g(1) - t46 * g(2)) * t45;
t227 = t49 * t221 - t226;
t229 = t211 * g(3);
t230 = t52 * t227 - t229;
t233 = cos(qJ(2,4));
t235 = t86 * t233 * t77;
t239 = -t70 * t67 + t71 * t68;
t240 = t233 * t239;
t243 = t68 * g(1) + g(2) * t70;
t248 = (t70 * g(1) - t68 * g(2)) * t67;
t249 = t71 * t243 - t248;
t251 = t233 * g(3);
t252 = t74 * t249 - t251;
t255 = cos(qJ(2,5));
t257 = t108 * t255 * t99;
t261 = -t92 * t89 + t93 * t90;
t262 = t255 * t261;
t265 = t90 * g(1) + g(2) * t92;
t270 = (t92 * g(1) - t90 * g(2)) * t89;
t271 = t93 * t265 - t270;
t273 = t255 * g(3);
t274 = t96 * t271 - t273;
t277 = cos(qJ(2,6));
t279 = t130 * t277 * t121;
t283 = -t114 * t111 + t115 * t112;
t284 = t277 * t283;
t287 = t112 * g(1) + g(2) * t114;
t292 = (t114 * g(1) - t112 * g(2)) * t111;
t293 = t115 * t287 - t292;
t295 = t277 * g(3);
t296 = t118 * t293 - t295;
t299 = t186 * t11 * t174 + t296 * t121 * t284 + t208 * t33 * t196 + t230 * t55 * t218 + t252 * t77 * t240 + t274 * t99 * t262 + t169 * t10 + t279 * t120 + t191 * t32 + t213 * t54 + t235 * t76 + t257 * t98;
t302 = t20 * t11 * t7;
t304 = t8 * g(3);
t305 = t167 * t183 + t304;
t309 = t42 * t33 * t29;
t311 = t30 * g(3);
t312 = t189 * t205 + t311;
t316 = t64 * t55 * t51;
t318 = t52 * g(3);
t319 = t211 * t227 + t318;
t323 = t86 * t77 * t73;
t325 = t74 * g(3);
t326 = t233 * t249 + t325;
t330 = t108 * t99 * t95;
t332 = t96 * g(3);
t333 = t255 * t271 + t332;
t337 = t130 * t121 * t117;
t339 = t118 * g(3);
t340 = t277 * t293 + t339;
t343 = t305 * t11 * t174 + t340 * t121 * t284 + t312 * t33 * t196 + t319 * t55 * t218 + t326 * t77 * t240 + t333 * t99 * t262 - t302 - t309 - t316 - t323 - t330 - t337;
t347 = t4 * t14 - t135;
t348 = t347 * t11;
t351 = t26 * t36 - t140;
t352 = t351 * t33;
t355 = t48 * t58 - t145;
t356 = t355 * t55;
t359 = t70 * t80 - t150;
t360 = t359 * t77;
t363 = t92 * t102 - t155;
t364 = t363 * t99;
t367 = t114 * t124 - t160;
t368 = t367 * t121;
t373 = -t5 * t177 + t182;
t375 = t167 * t373 - t304;
t379 = -t27 * t199 + t204;
t381 = t189 * t379 - t311;
t385 = -t49 * t221 + t226;
t387 = t211 * t385 - t318;
t391 = -t71 * t243 + t248;
t393 = t233 * t391 - t325;
t397 = -t93 * t265 + t270;
t399 = t255 * t397 - t332;
t403 = -t115 * t287 + t292;
t405 = t277 * t403 - t339;
t408 = t375 * t11 * t174 + t405 * t121 * t284 + t381 * t33 * t196 + t387 * t55 * t218 + t393 * t77 * t240 + t399 * t99 * t262 + t302 + t309 + t316 + t323 + t330 + t337;
t414 = t8 * t373 + t185;
t420 = t30 * t379 + t207;
t426 = t52 * t385 + t229;
t432 = t74 * t391 + t251;
t438 = t96 * t397 + t273;
t444 = t118 * t403 + t295;
t446 = t20 * t7 - t305 * t174 + t414 * t8 * t173 + t42 * t29 - t312 * t196 + t420 * t30 * t195 + t64 * t51 - t319 * t218 + t426 * t52 * t217 + t86 * t73 - t326 * t240 + t432 * t74 * t239 + t108 * t95 - t333 * t262 + t438 * t96 * t261 + t130 * t117 - t340 * t284 + t444 * t118 * t283;
t448 = cos(xP(5));
t449 = cos(xP(6));
t450 = t449 * t448;
t452 = sin(xP(4));
t453 = sin(xP(5));
t454 = t453 * t452;
t456 = cos(xP(4));
t457 = sin(xP(6));
t459 = t449 * t454 + t457 * t456;
t461 = t453 * t456;
t464 = -t449 * t461 + t457 * t452;
t466 = -g(1) * t450 - g(2) * t459 - g(3) * t464;
t468 = t457 * t448;
t472 = t449 * t456 - t457 * t454;
t476 = t449 * t452 + t457 * t461;
t478 = g(1) * t468 - g(2) * t472 - g(3) * t476;
t481 = t452 * t448;
t483 = t456 * t448;
t485 = -g(1) * t453 + g(2) * t481 - g(3) * t483;
t490 = t9 * t173;
t492 = t31 * t195;
t494 = t53 * t217;
t496 = t75 * t239;
t498 = t97 * t261;
t500 = t119 * t283;
t513 = -t11 * t7;
t517 = -t33 * t29;
t521 = -t55 * t51;
t525 = -t77 * t73;
t529 = -t99 * t95;
t533 = -t121 * t117;
t536 = t186 * t167 * t513 + t208 * t189 * t517 + t230 * t211 * t521 + t252 * t233 * t525 + t274 * t255 * t529 + t296 * t277 * t533 + t169 * t490 + t191 * t492 + t213 * t494 + t235 * t496 + t257 * t498 + t279 * t500;
t538 = t11 * t173;
t539 = t20 * t538;
t542 = t33 * t195;
t543 = t42 * t542;
t546 = t55 * t217;
t547 = t64 * t546;
t550 = t77 * t239;
t551 = t86 * t550;
t554 = t99 * t261;
t555 = t108 * t554;
t558 = t121 * t283;
t559 = t130 * t558;
t562 = t305 * t167 * t513 + t312 * t189 * t517 + t319 * t211 * t521 + t326 * t233 * t525 + t333 * t255 * t529 + t340 * t277 * t533 - t539 - t543 - t547 - t551 - t555 - t559;
t585 = t375 * t167 * t513 + t381 * t189 * t517 + t387 * t211 * t521 + t393 * t233 * t525 + t399 * t255 * t529 + t405 * t277 * t533 + t539 + t543 + t547 + t551 + t555 + t559;
t588 = -t167 * t7;
t593 = -t189 * t29;
t598 = -t211 * t51;
t603 = -t233 * t73;
t608 = -t255 * t95;
t613 = -t277 * t117;
t617 = t20 * t173 - t305 * t588 - t414 * t8 * t7 + t42 * t195 - t312 * t593 - t420 * t29 * t30 + t64 * t217 - t319 * t598 - t426 * t51 * t52 + t86 * t239 - t326 * t603 - t432 * t73 * t74 + t108 * t261 - t333 * t608 - t438 * t95 * t96 + t130 * t283 - t340 * t613 - t444 * t117 * t118;
t625 = t8 * t11;
t627 = t30 * t33;
t629 = t52 * t55;
t631 = t74 * t77;
t633 = t96 * t99;
t635 = t118 * t121;
t637 = t186 * t625 + t208 * t627 + t230 * t629 + t252 * t631 + t274 * t633 + t296 * t635;
t668 = -t340 * t118 - t414 * t167 - t420 * t189 - t426 * t211 - t432 * t233 - t438 * t255 - t444 * t277 - t312 * t30 - t305 * t8 - t319 * t52 - t326 * t74 - t333 * t96;
t676 = t453 * t449;
t678 = t457 * t453;
t681 = koppelP(1,3) * t448 - koppelP(1,1) * t676 + koppelP(1,2) * t678;
t685 = koppelP(1,2) * t449 + koppelP(1,1) * t457;
t687 = t685 * t452 + t456 * t681;
t688 = t687 * t173;
t689 = t11 * t9;
t690 = t20 * t689;
t695 = koppelP(2,3) * t448 - koppelP(2,1) * t676 + koppelP(2,2) * t678;
t699 = koppelP(2,2) * t449 + koppelP(2,1) * t457;
t701 = t699 * t452 + t456 * t695;
t702 = t701 * t195;
t703 = t33 * t31;
t704 = t42 * t703;
t709 = koppelP(3,3) * t448 - koppelP(3,1) * t676 + koppelP(3,2) * t678;
t713 = koppelP(3,2) * t449 + koppelP(3,1) * t457;
t715 = t713 * t452 + t456 * t709;
t716 = t715 * t217;
t717 = t55 * t53;
t718 = t64 * t717;
t723 = koppelP(4,3) * t448 - koppelP(4,1) * t676 + koppelP(4,2) * t678;
t727 = koppelP(4,2) * t449 + koppelP(4,1) * t457;
t729 = t727 * t452 + t456 * t723;
t730 = t729 * t239;
t731 = t77 * t75;
t732 = t86 * t731;
t737 = koppelP(5,3) * t448 - koppelP(5,1) * t676 + koppelP(5,2) * t678;
t741 = koppelP(5,2) * t449 + koppelP(5,1) * t457;
t743 = t741 * t452 + t456 * t737;
t744 = t743 * t261;
t745 = t99 * t97;
t746 = t108 * t745;
t751 = koppelP(6,3) * t448 - koppelP(6,1) * t676 + koppelP(6,2) * t678;
t755 = t449 * koppelP(6,2) + koppelP(6,1) * t457;
t757 = t755 * t452 + t456 * t751;
t758 = t283 * t757;
t759 = t121 * t119;
t760 = t130 * t759;
t764 = t137 * t689;
t766 = t142 * t703;
t768 = t147 * t717;
t770 = t152 * t731;
t772 = t157 * t745;
t774 = t162 * t759;
t780 = -t7 * t687;
t784 = -t681 * t452 + t685 * t456;
t786 = -t167 * t780 + t8 * t784;
t787 = t11 * t786;
t791 = -t29 * t701;
t795 = -t695 * t452 + t699 * t456;
t797 = -t189 * t791 + t30 * t795;
t798 = t33 * t797;
t802 = -t51 * t715;
t806 = -t709 * t452 + t713 * t456;
t808 = -t211 * t802 + t52 * t806;
t809 = t55 * t808;
t813 = -t73 * t729;
t817 = -t723 * t452 + t727 * t456;
t819 = -t233 * t813 + t74 * t817;
t820 = t77 * t819;
t824 = -t743 * t95;
t828 = -t737 * t452 + t741 * t456;
t830 = -t255 * t824 + t96 * t828;
t831 = t99 * t830;
t835 = -t757 * t117;
t839 = -t751 * t452 + t755 * t456;
t841 = t118 * t839 - t277 * t835;
t842 = t121 * t841;
t844 = -t279 * t119 * t758 - t169 * t9 * t688 - t191 * t31 * t702 - t213 * t53 * t716 - t235 * t75 * t730 - t257 * t97 * t744 + t186 * t787 + t208 * t798 + t230 * t809 + t252 * t820 + t274 * t831 + t296 * t842;
t846 = t21 * t688;
t848 = t43 * t702;
t850 = t65 * t716;
t852 = t87 * t730;
t854 = t109 * t744;
t856 = t131 * t758;
t858 = t305 * t787 + t312 * t798 + t319 * t809 + t326 * t820 + t333 * t831 + t340 * t842 + t846 + t848 + t850 + t852 + t854 + t856;
t861 = t347 * t689;
t863 = t351 * t703;
t865 = t355 * t717;
t867 = t359 * t731;
t869 = t363 * t745;
t871 = t367 * t759;
t881 = t375 * t787 + t381 * t798 + t387 * t809 + t393 * t820 + t399 * t831 + t405 * t842 - t846 - t848 - t850 - t852 - t854 - t856;
t919 = -t20 * t688 - t305 * t786 + t414 * (-t784 * t167 - t8 * t780) - t42 * t702 - t312 * t797 + t420 * (-t795 * t189 - t30 * t791) - t64 * t716 - t319 * t808 + t426 * (-t806 * t211 - t52 * t802) - t86 * t730 - t326 * t819 + t432 * (-t817 * t233 - t74 * t813) - t108 * t744 - t333 * t830 + t438 * (-t828 * t255 - t96 * t824) - t130 * t758 - t340 * t841 + t444 * (-t118 * t835 - t839 * t277);
t933 = MDP(2) * (-t690 * t688 - t704 * t702 - t718 * t716 - t732 * t730 - t746 * t744 - t760 * t758) + MDP(3) * (-t764 * t688 - t766 * t702 - t768 * t716 - t770 * t730 - t772 * t744 - t774 * t758) + MDP(9) * t844 + MDP(10) * t858 + MDP(11) * t844 + MDP(12) * (-t861 * t688 - t863 * t702 - t865 * t716 - t867 * t730 - t869 * t744 - t871 * t758) + MDP(13) * t881 + MDP(14) * t919 + MDP(21) * (t478 * t453 + t485 * t468) + MDP(22) * (t485 * t450 - t466 * t453) + MDP(23) * (-t478 * t450 - t466 * t468);
t955 = koppelP(1,3) * t453;
t956 = koppelP(1,1) * t450 - koppelP(1,2) * t468 + t955;
t958 = t167 * t688 - t956 * t8;
t959 = t11 * t958;
t966 = koppelP(2,3) * t453;
t967 = koppelP(2,1) * t450 - koppelP(2,2) * t468 + t966;
t969 = t189 * t702 - t967 * t30;
t970 = t33 * t969;
t977 = koppelP(3,3) * t453;
t978 = koppelP(3,1) * t450 - koppelP(3,2) * t468 + t977;
t980 = t211 * t716 - t978 * t52;
t981 = t55 * t980;
t988 = koppelP(4,3) * t453;
t989 = koppelP(4,1) * t450 - koppelP(4,2) * t468 + t988;
t991 = t233 * t730 - t989 * t74;
t992 = t77 * t991;
t999 = koppelP(5,3) * t453;
t1000 = koppelP(5,1) * t450 - koppelP(5,2) * t468 + t999;
t1002 = -t1000 * t96 + t255 * t744;
t1003 = t99 * t1002;
t1010 = koppelP(6,3) * t453;
t1011 = koppelP(6,1) * t450 - koppelP(6,2) * t468 + t1010;
t1013 = -t1011 * t118 + t277 * t758;
t1014 = t121 * t1013;
t1016 = -t279 * t119 * t835 - t169 * t9 * t780 - t191 * t31 * t791 - t213 * t53 * t802 - t235 * t75 * t813 - t257 * t97 * t824 + t274 * t1003 + t296 * t1014 + t186 * t959 + t208 * t970 + t230 * t981 + t252 * t992;
t1018 = t21 * t780;
t1020 = t43 * t791;
t1022 = t65 * t802;
t1024 = t87 * t813;
t1026 = t109 * t824;
t1028 = t131 * t835;
t1030 = t333 * t1003 + t340 * t1014 + t305 * t959 + t312 * t970 + t319 * t981 + t326 * t992 + t1018 + t1020 + t1022 + t1024 + t1026 + t1028;
t1047 = t399 * t1003 + t405 * t1014 + t375 * t959 + t381 * t970 + t387 * t981 + t393 * t992 - t1018 - t1020 - t1022 - t1024 - t1026 - t1028;
t1091 = -t20 * t780 - t305 * t958 + t414 * (t8 * t687 * t173 + t956 * t167) - t42 * t791 - t312 * t969 + t420 * (t30 * t701 * t195 + t967 * t189) - t64 * t802 - t319 * t980 + t426 * (t52 * t715 * t217 + t978 * t211) - t86 * t813 - t326 * t991 + t432 * (t74 * t729 * t239 + t989 * t233) - t108 * t824 - t333 * t1002 + t438 * (t96 * t743 * t261 + t1000 * t255) - t130 * t835 - t340 * t1013 + t444 * (t118 * t283 * t757 + t1011 * t277);
t1105 = MDP(2) * (-t690 * t780 - t704 * t791 - t718 * t802 - t732 * t813 - t746 * t824 - t760 * t835) + MDP(3) * (-t764 * t780 - t766 * t791 - t768 * t802 - t770 * t813 - t772 * t824 - t774 * t835) + MDP(9) * t1016 + MDP(10) * t1030 + MDP(11) * t1016 + MDP(12) * (-t861 * t780 - t863 * t791 - t865 * t802 - t867 * t813 - t869 * t824 - t871 * t835) + MDP(13) * t1047 + MDP(14) * t1091 + MDP(21) * (-t485 * t472 - t478 * t481) + MDP(22) * (t485 * t459 + t466 * t481) + MDP(23) * (-t478 * t459 + t466 * t472);
t1107 = t448 * koppelP(1,3) * t452;
t1110 = koppelP(1,1) * t454 + t456 * koppelP(1,2);
t1115 = (-koppelP(1,2) * t454 + koppelP(1,1) * t456) * t457;
t1116 = t449 * t1110 - t1107 + t1115;
t1121 = -t689 * t7 * t1116 + t689 * t173 * t956;
t1124 = t448 * koppelP(2,3) * t452;
t1127 = koppelP(2,1) * t454 + koppelP(2,2) * t456;
t1133 = (-t452 * koppelP(2,2) * t453 + koppelP(2,1) * t456) * t457;
t1134 = t449 * t1127 - t1124 + t1133;
t1139 = -t703 * t29 * t1134 + t703 * t195 * t967;
t1142 = t448 * koppelP(3,3) * t452;
t1145 = koppelP(3,1) * t454 + koppelP(3,2) * t456;
t1151 = (-t452 * koppelP(3,2) * t453 + koppelP(3,1) * t456) * t457;
t1152 = t449 * t1145 - t1142 + t1151;
t1157 = -t717 * t51 * t1152 + t717 * t217 * t978;
t1160 = t448 * koppelP(4,3) * t452;
t1163 = koppelP(4,1) * t454 + koppelP(4,2) * t456;
t1169 = (-t452 * koppelP(4,2) * t453 + koppelP(4,1) * t456) * t457;
t1170 = t449 * t1163 - t1160 + t1169;
t1175 = -t731 * t73 * t1170 + t731 * t239 * t989;
t1178 = t448 * koppelP(5,3) * t452;
t1181 = koppelP(5,1) * t454 + koppelP(5,2) * t456;
t1187 = (-t452 * koppelP(5,2) * t453 + koppelP(5,1) * t456) * t457;
t1188 = t449 * t1181 - t1178 + t1187;
t1193 = t745 * t261 * t1000 - t745 * t95 * t1188;
t1196 = t448 * koppelP(6,3) * t452;
t1199 = koppelP(6,1) * t454 + koppelP(6,2) * t456;
t1205 = (-t452 * koppelP(6,2) * t453 + koppelP(6,1) * t456) * t457;
t1206 = t449 * t1199 - t1196 + t1205;
t1211 = t759 * t283 * t1011 - t759 * t117 * t1206;
t1229 = t588 * t11 * t956 - t538 * t1116 * t167;
t1237 = -t542 * t1134 * t189 + t593 * t33 * t967;
t1245 = -t546 * t1152 * t211 + t598 * t55 * t978;
t1253 = -t550 * t1170 * t233 + t603 * t77 * t989;
t1261 = t608 * t99 * t1000 - t554 * t1188 * t255;
t1269 = t613 * t121 * t1011 - t558 * t1206 * t277;
t1271 = t108 * t255 * t1193 + t20 * t167 * t1121 + t42 * t189 * t1139 + t64 * t211 * t1157 + t86 * t233 * t1175 + t130 * t277 * t1211 + t186 * t1229 + t208 * t1237 + t230 * t1245 + t252 * t1253 + t274 * t1261 + t296 * t1269;
t1274 = t20 * t8 * t1121;
t1277 = t42 * t30 * t1139;
t1280 = t64 * t52 * t1157;
t1283 = t86 * t74 * t1175;
t1286 = t108 * t96 * t1193;
t1289 = t130 * t118 * t1211;
t1291 = t305 * t1229 + t312 * t1237 + t319 * t1245 + t326 * t1253 + t333 * t1261 + t340 * t1269 - t1274 - t1277 - t1280 - t1283 - t1286 - t1289;
t1308 = t375 * t1229 + t381 * t1237 + t387 * t1245 + t393 * t1253 + t399 * t1261 + t405 * t1269 + t1274 + t1277 + t1280 + t1283 + t1286 + t1289;
t1322 = t448 * (t449 * koppelP(1,1) - koppelP(1,2) * t457) + t955;
t1345 = t448 * (t449 * koppelP(2,1) - koppelP(2,2) * t457) + t966;
t1368 = t448 * (t449 * koppelP(3,1) - koppelP(3,2) * t457) + t977;
t1391 = t448 * (t449 * koppelP(4,1) - koppelP(4,2) * t457) + t988;
t1414 = t448 * (t449 * koppelP(5,1) - koppelP(5,2) * t457) + t999;
t1437 = t448 * (t449 * koppelP(6,1) - koppelP(6,2) * t457) + t1010;
t1448 = t20 * t8 * qJ(3,1) * t1121 - t305 * qJ(3,1) * t1229 + t414 * t8 * (t5 * (t2 * (-t449 * t1110 + t1107 - t1115) + t1322 * t4) + t1 * (t1116 * t4 + t2 * t1322)) + t42 * t30 * qJ(3,2) * t1139 - t312 * qJ(3,2) * t1237 + t420 * t30 * (t27 * (t24 * (-t449 * t1127 + t1124 - t1133) + t1345 * t26) + t23 * (t1134 * t26 + t24 * t1345)) + t64 * t52 * qJ(3,3) * t1157 - t319 * qJ(3,3) * t1245 + t426 * t52 * (t49 * (t46 * (-t449 * t1145 + t1142 - t1151) + t1368 * t48) + t45 * (t1152 * t48 + t46 * t1368)) + t86 * t74 * qJ(3,4) * t1175 - t326 * qJ(3,4) * t1253 + t432 * t74 * (t71 * (t68 * (-t449 * t1163 + t1160 - t1169) + t1391 * t70) + t67 * (t1170 * t70 + t68 * t1391)) + t108 * t96 * qJ(3,5) * t1193 - t333 * qJ(3,5) * t1261 + t438 * t96 * (t93 * (t90 * (-t449 * t1181 + t1178 - t1187) + t1414 * t92) + t89 * (t1188 * t92 + t90 * t1414)) + t130 * t118 * qJ(3,6) * t1211 - t340 * qJ(3,6) * t1269 + t444 * t118 * (t115 * (t112 * (-t449 * t1199 + t1196 - t1205) + t1437 * t114) + t111 * (t112 * t1437 + t1206 * t114));
t1462 = MDP(2) * (t108 * t1193 + t20 * t1121 + t42 * t1139 + t64 * t1157 + t86 * t1175 + t130 * t1211) + MDP(3) * (t137 * t1121 + t142 * t1139 + t147 * t1157 + t152 * t1175 + t157 * t1193 + t162 * t1211) + MDP(9) * t1271 + MDP(10) * t1291 + MDP(11) * t1271 + MDP(12) * (t347 * t1121 + t351 * t1139 + t355 * t1157 + t359 * t1175 + t363 * t1193 + t367 * t1211) + MDP(13) * t1308 + MDP(14) * t1448 + MDP(21) * (-t485 * t476 + t478 * t483) + MDP(22) * (t485 * t464 - t466 * t483) + MDP(23) * (-t478 * t464 + t466 * t476);
unknown(1,1) = MDP(2) * (t21 * t10 + t109 * t98 + t131 * t120 + t43 * t32 + t65 * t54 + t87 * t76) + MDP(3) * (t138 * t10 + t163 * t120 + t143 * t32 + t148 * t54 + t153 * t76 + t158 * t98) + MDP(9) * t299 + MDP(10) * t343 + MDP(11) * t299 + MDP(12) * (t348 * t10 + t368 * t120 + t352 * t32 + t356 * t54 + t360 * t76 + t364 * t98) + MDP(13) * t408 + MDP(14) * t446 + MDP(24) * (t466 * t450 + t485 * t453 - t478 * t468);
unknown(2,1) = MDP(2) * (t109 * t498 + t131 * t500 + t21 * t490 + t43 * t492 + t65 * t494 + t87 * t496) + MDP(3) * (t138 * t490 + t143 * t492 + t148 * t494 + t153 * t496 + t158 * t498 + t163 * t500) + MDP(9) * t536 + MDP(10) * t562 + MDP(11) * t536 + MDP(12) * (t348 * t490 + t352 * t492 + t356 * t494 + t360 * t496 + t364 * t498 + t368 * t500) + MDP(13) * t585 + MDP(14) * t617 + MDP(24) * (t466 * t459 + t478 * t472 - t485 * t481);
unknown(3,1) = MDP(9) * t637 + MDP(10) * (t305 * t625 + t312 * t627 + t319 * t629 + t326 * t631 + t333 * t633 + t340 * t635) + MDP(11) * t637 + MDP(13) * (t375 * t625 + t381 * t627 + t387 * t629 + t393 * t631 + t399 * t633 + t405 * t635) + MDP(14) * t668 + MDP(24) * (t466 * t464 + t478 * t476 + t485 * t483);
unknown(4,1) = t933;
unknown(5,1) = t1105;
unknown(6,1) = t1462;
taugX  = unknown;
