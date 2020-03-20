% Calculate Gravitation load for parallel robot
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
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 23:28
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1P4A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 1
% StartTime: 2020-03-12 23:26:55
% EndTime: 2020-03-12 23:26:56
% DurationCPUTime: 0.34s
% Computational Cost: add. (2141->361), mult. (4125->777), div. (126->12), fcn. (2783->42), ass. (0->271)
unknown=NaN(6,1);
t1 = sin(legFrame(1,3));
t2 = cos(qJ(1,1));
t4 = sin(qJ(1,1));
t5 = cos(legFrame(1,3));
t7 = -t2 * t1 - t5 * t4;
t8 = sin(qJ(2,1));
t9 = 0.1e1 / t8;
t11 = 0.1e1 / qJ(3,1);
t12 = rSges(3,3) + qJ(3,1);
t14 = m(2) * rSges(2,2);
t15 = m(3) * t12 - t14;
t18 = -g(1) * t1 + g(2) * t5;
t23 = rSges(2,1) * m(2) + m(3) * rSges(3,1);
t24 = t23 * t18;
t25 = cos(qJ(2,1));
t30 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - m(3) * rSges(3,2);
t33 = g(1) * t5 + g(2) * t1;
t41 = t23 * t33;
t48 = t2 * (-rSges(1,1) * t18 * m(1) - t8 * t18 * t15 - t25 * t24 + t33 * t30) + (rSges(1,1) * t33 * m(1) + t8 * t33 * t15 + t30 * t18 + t25 * t41) * t4;
t49 = t48 * t11;
t53 = -t4 * t1 + t5 * t2;
t56 = -m(3) * t12 + t14;
t61 = t23 * g(3);
t69 = t25 * (t4 * t18 * t56 + t2 * t33 * t56 - t61) + t8 * (g(3) * t56 + t2 * t41 + t4 * t24);
t78 = t8 * t4 * t18 + t8 * t2 * t33 - t25 * g(3);
t79 = t78 * m(3);
t81 = sin(legFrame(2,3));
t82 = cos(qJ(1,2));
t84 = sin(qJ(1,2));
t85 = cos(legFrame(2,3));
t87 = -t82 * t81 - t85 * t84;
t88 = sin(qJ(2,2));
t89 = 0.1e1 / t88;
t91 = 0.1e1 / qJ(3,2);
t92 = rSges(3,3) + qJ(3,2);
t94 = m(3) * t92 - t14;
t97 = -g(1) * t81 + g(2) * t85;
t100 = t23 * t97;
t101 = cos(qJ(2,2));
t105 = g(1) * t85 + g(2) * t81;
t113 = t23 * t105;
t120 = t82 * (-rSges(1,1) * t97 * m(1) - t88 * t97 * t94 - t101 * t100 + t105 * t30) + (rSges(1,1) * t105 * m(1) + t88 * t105 * t94 + t101 * t113 + t30 * t97) * t84;
t121 = t120 * t91;
t125 = -t84 * t81 + t85 * t82;
t128 = -m(3) * t92 + t14;
t140 = t101 * (t82 * t105 * t128 + t84 * t97 * t128 - t61) + t88 * (g(3) * t128 + t84 * t100 + t82 * t113);
t149 = t88 * t82 * t105 + t88 * t84 * t97 - t101 * g(3);
t150 = t149 * m(3);
t152 = sin(legFrame(3,3));
t153 = cos(qJ(1,3));
t155 = sin(qJ(1,3));
t156 = cos(legFrame(3,3));
t158 = -t153 * t152 - t156 * t155;
t159 = sin(qJ(2,3));
t160 = 0.1e1 / t159;
t162 = 0.1e1 / qJ(3,3);
t163 = rSges(3,3) + qJ(3,3);
t165 = m(3) * t163 - t14;
t168 = -g(1) * t152 + g(2) * t156;
t171 = t23 * t168;
t172 = cos(qJ(2,3));
t176 = g(1) * t156 + g(2) * t152;
t184 = t23 * t176;
t191 = t153 * (-rSges(1,1) * t168 * m(1) - t159 * t168 * t165 - t172 * t171 + t176 * t30) + (rSges(1,1) * t176 * m(1) + t159 * t176 * t165 + t30 * t168 + t172 * t184) * t155;
t192 = t191 * t162;
t196 = -t155 * t152 + t156 * t153;
t199 = -m(3) * t163 + t14;
t211 = t172 * (t153 * t176 * t199 + t155 * t168 * t199 - t61) + t159 * (g(3) * t199 + t153 * t184 + t155 * t171);
t220 = t159 * t153 * t176 + t159 * t155 * t168 - t172 * g(3);
t221 = t220 * m(3);
t223 = sin(legFrame(4,3));
t224 = cos(qJ(1,4));
t226 = sin(qJ(1,4));
t227 = cos(legFrame(4,3));
t229 = -t224 * t223 - t227 * t226;
t230 = sin(qJ(2,4));
t231 = 0.1e1 / t230;
t233 = 0.1e1 / qJ(3,4);
t234 = rSges(3,3) + qJ(3,4);
t236 = m(3) * t234 - t14;
t239 = -g(1) * t223 + g(2) * t227;
t242 = t23 * t239;
t243 = cos(qJ(2,4));
t247 = g(1) * t227 + g(2) * t223;
t255 = t23 * t247;
t262 = t224 * (-rSges(1,1) * t239 * m(1) - t230 * t239 * t236 - t243 * t242 + t247 * t30) + (rSges(1,1) * t247 * m(1) + t230 * t247 * t236 + t30 * t239 + t243 * t255) * t226;
t263 = t262 * t233;
t267 = -t226 * t223 + t227 * t224;
t270 = -m(3) * t234 + t14;
t282 = t243 * (t224 * t247 * t270 + t226 * t239 * t270 - t61) + t230 * (g(3) * t270 + t224 * t255 + t226 * t242);
t291 = t230 * t224 * t247 + t230 * t226 * t239 - t243 * g(3);
t292 = t291 * m(3);
t294 = sin(legFrame(5,3));
t295 = cos(qJ(1,5));
t297 = sin(qJ(1,5));
t298 = cos(legFrame(5,3));
t300 = -t295 * t294 - t298 * t297;
t301 = sin(qJ(2,5));
t302 = 0.1e1 / t301;
t304 = 0.1e1 / qJ(3,5);
t305 = rSges(3,3) + qJ(3,5);
t307 = m(3) * t305 - t14;
t310 = -g(1) * t294 + g(2) * t298;
t313 = t23 * t310;
t314 = cos(qJ(2,5));
t318 = g(1) * t298 + g(2) * t294;
t326 = t23 * t318;
t333 = t295 * (-rSges(1,1) * t310 * m(1) - t301 * t310 * t307 + t318 * t30 - t314 * t313) + (rSges(1,1) * t318 * m(1) + t301 * t318 * t307 + t30 * t310 + t314 * t326) * t297;
t334 = t333 * t304;
t338 = -t297 * t294 + t298 * t295;
t341 = -m(3) * t305 + t14;
t353 = t314 * (t295 * t318 * t341 + t297 * t310 * t341 - t61) + t301 * (g(3) * t341 + t295 * t326 + t297 * t313);
t362 = t301 * t295 * t318 + t301 * t297 * t310 - t314 * g(3);
t363 = t362 * m(3);
t365 = sin(legFrame(6,3));
t366 = cos(qJ(1,6));
t368 = sin(qJ(1,6));
t369 = cos(legFrame(6,3));
t371 = -t366 * t365 - t369 * t368;
t372 = sin(qJ(2,6));
t373 = 0.1e1 / t372;
t375 = 0.1e1 / qJ(3,6);
t376 = rSges(3,3) + qJ(3,6);
t378 = m(3) * t376 - t14;
t381 = -g(1) * t365 + g(2) * t369;
t384 = t23 * t381;
t385 = cos(qJ(2,6));
t389 = g(1) * t369 + g(2) * t365;
t397 = t23 * t389;
t404 = t366 * (-rSges(1,1) * t381 * m(1) - t372 * t381 * t378 + t389 * t30 - t385 * t384) + (rSges(1,1) * t389 * m(1) + t372 * t389 * t378 + t30 * t381 + t385 * t397) * t368;
t405 = t404 * t375;
t409 = -t368 * t365 + t369 * t366;
t412 = -m(3) * t376 + t14;
t424 = t385 * (t366 * t389 * t412 + t368 * t381 * t412 - t61) + t372 * (g(3) * t412 + t366 * t397 + t368 * t384);
t433 = t372 * t366 * t389 + t372 * t368 * t381 - t385 * g(3);
t434 = t433 * m(3);
t437 = t49 * t9 * t7 + t69 * t11 * t25 * t53 - t79 * t8 * t53 + t121 * t89 * t87 + t140 * t91 * t101 * t125 - t150 * t88 * t125 + t192 * t160 * t158 + t211 * t162 * t172 * t196 - t221 * t159 * t196 + t263 * t231 * t229 + t282 * t233 * t243 * t267 - t292 * t230 * t267 + t334 * t302 * t300 + t353 * t304 * t314 * t338 - t363 * t301 * t338 + t405 * t373 * t371 + t424 * t375 * t385 * t409 - t434 * t372 * t409 - m(4) * g(1);
t481 = t49 * t9 * t53 - t69 * t25 * t11 * t7 + t79 * t7 * t8 + t121 * t89 * t125 - t140 * t101 * t87 * t91 + t150 * t87 * t88 + t192 * t160 * t196 - t211 * t172 * t158 * t162 + t221 * t158 * t159 + t263 * t231 * t267 - t282 * t243 * t229 * t233 + t292 * t229 * t230 + t334 * t302 * t338 - t353 * t314 * t300 * t304 + t363 * t300 * t301 + t405 * t373 * t409 - t424 * t385 * t371 * t375 + t434 * t371 * t372 - m(4) * g(2);
t507 = t69 * t8 * t11 + t78 * m(3) * t25 + t140 * t88 * t91 + t149 * m(3) * t101 + t211 * t159 * t162 + t220 * m(3) * t172 + t282 * t230 * t233 + t291 * m(3) * t243 + t353 * t301 * t304 + t362 * m(3) * t314 + t424 * t372 * t375 + t433 * m(3) * t385 - m(4) * g(3);
t508 = sin(xP(6));
t510 = cos(xP(6));
t512 = -koppelP(1,2) * t508 + t510 * koppelP(1,1);
t513 = sin(xP(5));
t515 = cos(xP(5));
t516 = koppelP(1,3) * t515;
t517 = t513 * t512 - t516;
t518 = cos(xP(4));
t520 = sin(xP(4));
t523 = koppelP(1,1) * t508 + koppelP(1,2) * t510;
t524 = t523 * t520;
t525 = t518 * t517 - t524;
t527 = t11 * t9;
t528 = t48 * t527;
t531 = -t25 * t7;
t534 = t523 * t518;
t535 = t520 * t517 + t534;
t548 = -koppelP(2,2) * t508 + t510 * koppelP(2,1);
t550 = koppelP(2,3) * t515;
t551 = t513 * t548 - t550;
t555 = koppelP(2,1) * t508 + koppelP(2,2) * t510;
t556 = t555 * t520;
t557 = t518 * t551 - t556;
t559 = t91 * t89;
t560 = t120 * t559;
t563 = -t101 * t87;
t566 = t555 * t518;
t567 = t520 * t551 + t566;
t580 = -koppelP(3,2) * t508 + t510 * koppelP(3,1);
t582 = koppelP(3,3) * t515;
t583 = t513 * t580 - t582;
t587 = koppelP(3,1) * t508 + koppelP(3,2) * t510;
t588 = t587 * t520;
t589 = t518 * t583 - t588;
t591 = t162 * t160;
t592 = t191 * t591;
t595 = -t172 * t158;
t598 = t587 * t518;
t599 = t520 * t583 + t598;
t612 = -koppelP(4,2) * t508 + t510 * koppelP(4,1);
t614 = koppelP(4,3) * t515;
t615 = t513 * t612 - t614;
t619 = koppelP(4,1) * t508 + koppelP(4,2) * t510;
t620 = t619 * t520;
t621 = t518 * t615 - t620;
t623 = t233 * t231;
t624 = t262 * t623;
t627 = -t243 * t229;
t630 = t619 * t518;
t631 = t520 * t615 + t630;
t644 = -koppelP(5,2) * t508 + t510 * koppelP(5,1);
t646 = koppelP(5,3) * t515;
t647 = t513 * t644 - t646;
t651 = koppelP(5,1) * t508 + koppelP(5,2) * t510;
t652 = t651 * t520;
t653 = t518 * t647 - t652;
t655 = t304 * t302;
t656 = t333 * t655;
t659 = -t314 * t300;
t662 = t651 * t518;
t663 = t520 * t647 + t662;
t676 = -koppelP(6,2) * t508 + t510 * koppelP(6,1);
t678 = koppelP(6,3) * t515;
t679 = t513 * t676 - t678;
t683 = koppelP(6,1) * t508 + t510 * koppelP(6,2);
t684 = t683 * t520;
t685 = t518 * t679 - t684;
t687 = t375 * t373;
t688 = t404 * t687;
t691 = -t385 * t371;
t694 = t683 * t518;
t695 = t520 * t679 + t694;
t706 = t513 * g(2);
t708 = g(3) * rSges(4,2);
t711 = rSges(4,2) * t508;
t719 = g(3) * t513;
t726 = t508 * g(2);
t734 = t528 * t53 * t525 + t69 * (t531 * t11 * t525 + t8 * t11 * t535) - t78 * m(3) * (-t7 * t8 * t525 - t25 * t535) + t560 * t125 * t557 + t140 * (t563 * t91 * t557 + t88 * t91 * t567) - t149 * m(3) * (-t87 * t88 * t557 - t101 * t567) + t592 * t196 * t589 + t211 * (t159 * t162 * t599 + t595 * t162 * t589) - t220 * m(3) * (-t158 * t159 * t589 - t172 * t599) + t624 * t267 * t621 + t282 * (t230 * t233 * t631 + t627 * t233 * t621) - t291 * m(3) * (-t229 * t230 * t621 - t243 * t631) + t656 * t338 * t653 + t353 * (t301 * t304 * t663 + t659 * t304 * t653) - t362 * m(3) * (-t300 * t301 * t653 - t314 * t663) + t688 * t409 * t685 + t424 * (t372 * t375 * t695 + t691 * t375 * t685) - t433 * m(3) * (-t371 * t372 * t685 - t385 * t695) + (t518 * (t510 * (-rSges(4,1) * t706 - t708) + t711 * t706 + rSges(4,3) * t515 * g(2) - t508 * rSges(4,1) * g(3)) + t520 * (t510 * (g(2) * rSges(4,2) - rSges(4,1) * t719) + t708 * t508 * t513 + rSges(4,1) * t726 + rSges(4,3) * g(3) * t515)) * m(4);
t736 = -t513 * t512 + t516;
t738 = t518 * t736 + t524;
t741 = t53 * t738;
t742 = t11 * t25;
t744 = t510 * t515;
t746 = t508 * t515;
t749 = -koppelP(1,3) * t513 - koppelP(1,1) * t744 + koppelP(1,2) * t746;
t760 = -t513 * t548 + t550;
t762 = t518 * t760 + t556;
t765 = t125 * t762;
t766 = t91 * t101;
t771 = -koppelP(2,3) * t513 - koppelP(2,1) * t744 + koppelP(2,2) * t746;
t782 = -t513 * t580 + t582;
t784 = t518 * t782 + t588;
t787 = t196 * t784;
t788 = t162 * t172;
t793 = -koppelP(3,3) * t513 - koppelP(3,1) * t744 + koppelP(3,2) * t746;
t804 = -t513 * t612 + t614;
t806 = t518 * t804 + t620;
t809 = t267 * t806;
t810 = t233 * t243;
t815 = -koppelP(4,3) * t513 - koppelP(4,1) * t744 + koppelP(4,2) * t746;
t826 = -t513 * t644 + t646;
t828 = t518 * t826 + t652;
t831 = t338 * t828;
t832 = t304 * t314;
t837 = -koppelP(5,3) * t513 - koppelP(5,1) * t744 + koppelP(5,2) * t746;
t848 = -t513 * t676 + t678;
t850 = t518 * t848 + t684;
t853 = t409 * t850;
t854 = t375 * t385;
t859 = -koppelP(6,3) * t513 - koppelP(6,1) * t744 + koppelP(6,2) * t746;
t871 = rSges(4,1) * t510 - t711;
t882 = g(1) * rSges(4,1);
t887 = t528 * t7 * t738 + t69 * (t8 * t11 * t749 + t742 * t741) - t78 * m(3) * (-t25 * t749 + t8 * t741) + t560 * t87 * t762 + t140 * (t88 * t91 * t771 + t766 * t765) - t149 * m(3) * (-t101 * t771 + t88 * t765) + t592 * t158 * t784 + t211 * (t159 * t162 * t793 + t788 * t787) - t220 * m(3) * (t159 * t787 - t172 * t793) + t624 * t229 * t806 + t282 * (t230 * t233 * t815 + t810 * t809) - t291 * m(3) * (t230 * t809 - t243 * t815) + t656 * t300 * t828 + t353 * (t301 * t304 * t837 + t832 * t831) - t362 * m(3) * (t301 * t831 - t314 * t837) + t688 * t371 * t850 + t424 * (t372 * t375 * t859 + t854 * t853) - t433 * m(3) * (t372 * t853 - t385 * t859) - m(4) * (t518 * g(1) * (-t871 * t513 + rSges(4,3) * t515) - t515 * t871 * g(3) + g(1) * rSges(4,2) * t520 * t510 + t882 * t508 * t520 - rSges(4,3) * t719);
t889 = t520 * t736 - t534;
t896 = t53 * t889;
t909 = t520 * t760 - t566;
t916 = t125 * t909;
t929 = t520 * t782 - t598;
t936 = t196 * t929;
t949 = t520 * t804 - t630;
t956 = t267 * t949;
t969 = t520 * t826 - t662;
t976 = t338 * t969;
t989 = t520 * t848 - t694;
t996 = t409 * t989;
t1030 = t48 * (-t527 * t53 * t749 + t527 * t7 * t889) + t69 * (-t531 * t11 * t749 + t742 * t896) - t78 * m(3) * (t7 * t8 * t749 + t8 * t896) + t120 * (-t559 * t125 * t771 + t559 * t87 * t909) + t140 * (-t563 * t91 * t771 + t766 * t916) - t149 * m(3) * (t87 * t88 * t771 + t88 * t916) + t191 * (t591 * t158 * t929 - t591 * t196 * t793) + t211 * (-t595 * t162 * t793 + t788 * t936) - t220 * m(3) * (t158 * t159 * t793 + t159 * t936) + t262 * (t623 * t229 * t949 - t623 * t267 * t815) + t282 * (-t627 * t233 * t815 + t810 * t956) - t291 * m(3) * (t229 * t230 * t815 + t230 * t956) + t333 * (t655 * t300 * t969 - t655 * t338 * t837) + t353 * (-t659 * t304 * t837 + t832 * t976) - t362 * m(3) * (t300 * t301 * t837 + t301 * t976) + t404 * (t687 * t371 * t989 - t687 * t409 * t859) + t424 * (-t691 * t375 * t859 + t854 * t996) - t433 * m(3) * (t371 * t372 * t859 + t372 * t996) + (t515 * (-rSges(4,3) * t520 * g(1) - t510 * rSges(4,1) * g(2) + rSges(4,2) * t726) + t510 * (t520 * rSges(4,1) * t513 + rSges(4,2) * t518) * g(1) - rSges(4,2) * g(1) * t508 * t520 * t513 + t882 * t518 * t508 - rSges(4,3) * t706) * m(4);
unknown(1,1) = t437;
unknown(2,1) = t481;
unknown(3,1) = t507;
unknown(4,1) = t734;
unknown(5,1) = t887;
unknown(6,1) = t1030;
taugX  = unknown;
