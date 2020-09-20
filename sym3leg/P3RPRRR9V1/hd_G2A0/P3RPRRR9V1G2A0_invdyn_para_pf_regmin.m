% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR9V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x15]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR9V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:50
% EndTime: 2020-08-06 18:52:06
% DurationCPUTime: 15.54s
% Computational Cost: add. (83836->735), mult. (128112->1230), div. (8889->17), fcn. (97065->87), ass. (0->570)
t705 = -2 * pkin(1);
t703 = 2 * pkin(3);
t355 = cos(pkin(7));
t314 = t355 * pkin(2);
t280 = t314 + pkin(1);
t356 = pkin(5) + qJ(2,3);
t340 = pkin(6) + t356;
t317 = 0.1e1 / t340;
t318 = 0.1e1 / t340 ^ 2;
t319 = t317 * t318;
t357 = pkin(5) + qJ(2,2);
t341 = pkin(6) + t357;
t320 = 0.1e1 / t341;
t321 = 0.1e1 / t341 ^ 2;
t322 = t320 * t321;
t358 = pkin(5) + qJ(2,1);
t342 = pkin(6) + t358;
t323 = 0.1e1 / t342;
t324 = 0.1e1 / t342 ^ 2;
t325 = t323 * t324;
t347 = pkin(7) + qJ(3,1);
t305 = cos(t347);
t663 = pkin(3) * t305;
t174 = t280 + t663;
t346 = pkin(7) + qJ(3,2);
t304 = cos(t346);
t664 = pkin(3) * t304;
t173 = t280 + t664;
t345 = pkin(7) + qJ(3,3);
t303 = cos(t345);
t665 = pkin(3) * t303;
t172 = t280 + t665;
t361 = legFrame(1,2);
t328 = sin(t361);
t331 = cos(t361);
t177 = t331 * g(1) - t328 * g(2);
t370 = sin(qJ(1,1));
t376 = cos(qJ(1,1));
t156 = g(3) * t376 + t177 * t370;
t360 = legFrame(2,2);
t327 = sin(t360);
t330 = cos(t360);
t176 = t330 * g(1) - t327 * g(2);
t368 = sin(qJ(1,2));
t374 = cos(qJ(1,2));
t155 = g(3) * t374 + t176 * t368;
t359 = legFrame(3,2);
t326 = sin(t359);
t329 = cos(t359);
t175 = t329 * g(1) - t326 * g(2);
t366 = sin(qJ(1,3));
t372 = cos(qJ(1,3));
t154 = g(3) * t372 + t175 * t366;
t704 = 2 * pkin(1);
t287 = 0.1e1 / t303;
t290 = 0.1e1 / t304;
t293 = 0.1e1 / t305;
t288 = 0.1e1 / t303 ^ 2;
t291 = 0.1e1 / t304 ^ 2;
t294 = 0.1e1 / t305 ^ 2;
t343 = t355 ^ 2;
t702 = 0.2e1 * t343;
t701 = 0.4e1 * t343;
t371 = cos(qJ(3,3));
t350 = t371 ^ 2;
t700 = 0.2e1 * t350;
t373 = cos(qJ(3,2));
t351 = t373 ^ 2;
t699 = 0.2e1 * t351;
t375 = cos(qJ(3,1));
t352 = t375 ^ 2;
t698 = 0.2e1 * t352;
t697 = -4 * pkin(5) - 4 * pkin(6);
t696 = -g(1) / 0.4e1;
t695 = g(1) / 0.4e1;
t694 = -g(2) / 0.4e1;
t693 = g(2) / 0.4e1;
t692 = g(3) / 0.2e1;
t289 = t287 * t288;
t297 = sin(t345);
t363 = xDDP(2);
t364 = xDDP(1);
t393 = 1 / pkin(3);
t379 = xDP(2);
t380 = xDP(1);
t160 = t326 * t380 + t329 * t379;
t157 = t160 ^ 2;
t394 = 1 / pkin(3) ^ 2;
t611 = t157 * t394;
t100 = t289 * t297 * t611 + (t326 * t364 + t329 * t363) * t393 * t287;
t378 = xDP(3);
t444 = t326 * t379 - t329 * t380;
t136 = -t366 * t444 + t372 * t378;
t583 = t297 * t160;
t103 = t136 * t303 + t583;
t499 = t160 * t288 * t317;
t462 = t393 * t499;
t453 = t103 * t462;
t691 = -t100 * t356 + t453 * t705;
t292 = t290 * t291;
t298 = sin(t346);
t161 = t327 * t380 + t330 * t379;
t158 = t161 ^ 2;
t610 = t158 * t394;
t101 = t298 * t292 * t610 + (t327 * t364 + t330 * t363) * t393 * t290;
t443 = t327 * t379 - t330 * t380;
t137 = -t368 * t443 + t374 * t378;
t582 = t298 * t161;
t104 = t137 * t304 + t582;
t498 = t161 * t291 * t320;
t461 = t393 * t498;
t452 = t104 * t461;
t690 = -t101 * t357 + t452 * t705;
t295 = t293 * t294;
t299 = sin(t347);
t162 = t328 * t380 + t331 * t379;
t159 = t162 ^ 2;
t609 = t159 * t394;
t102 = t299 * t295 * t609 + (t328 * t364 + t331 * t363) * t393 * t293;
t442 = t328 * t379 - t331 * t380;
t138 = -t370 * t442 + t376 * t378;
t581 = t299 * t162;
t105 = t138 * t305 + t581;
t497 = t162 * t294 * t323;
t460 = t393 * t497;
t451 = t105 * t460;
t689 = -t102 * t358 + t451 * t705;
t354 = sin(pkin(7));
t365 = sin(qJ(3,3));
t570 = t354 * t365;
t167 = t355 * t371 - t570;
t164 = 0.1e1 / t167;
t338 = pkin(1) * t378;
t562 = (-pkin(1) * t444 + t340 * t378) * t366 + (t340 * t444 + t338) * t372;
t566 = t371 * t365;
t315 = pkin(1) * t354;
t567 = t371 * (-t365 * pkin(3) + t315);
t608 = t160 * t365;
t620 = t136 * t350;
t61 = ((t136 * t371 + t608) * pkin(2) + (0.2e1 * t160 * t566 - t136 + 0.2e1 * t620) * pkin(3)) * t343 + (t562 * t371 + pkin(1) * t608 + ((-t136 * t365 + t160 * t371) * pkin(2) + (-0.2e1 * t136 * t566 + t160 * t700 - t160) * pkin(3)) * t354) * t355 - pkin(3) * t620 + t160 * t567 - t562 * t570 + pkin(3) * t136;
t629 = t61 * t164;
t537 = t319 * t629;
t638 = (t629 + (-t280 * t287 - pkin(3)) * t103) * t319;
t428 = (-t537 - t638) * t287;
t577 = t326 * t366;
t142 = t329 * t297 - t303 * t577;
t592 = t287 * t317;
t514 = t142 * t592;
t106 = t363 * t514;
t574 = t329 * t366;
t143 = t326 * t297 + t303 * t574;
t513 = t143 * t592;
t107 = t364 * t513;
t133 = t157 * t393 * t289 * t317;
t362 = xDDP(3);
t580 = t317 * t372;
t205 = t362 * t580;
t456 = t106 + t107 + t133 + t205;
t49 = t103 * t428 + t456;
t688 = t49 * pkin(1);
t367 = sin(qJ(3,2));
t569 = t354 * t367;
t168 = t355 * t373 - t569;
t165 = 0.1e1 / t168;
t561 = (-pkin(1) * t443 + t341 * t378) * t368 + (t341 * t443 + t338) * t374;
t564 = t373 * t367;
t565 = t373 * (-t367 * pkin(3) + t315);
t607 = t161 * t367;
t619 = t137 * t351;
t62 = ((t137 * t373 + t607) * pkin(2) + (0.2e1 * t161 * t564 - t137 + 0.2e1 * t619) * pkin(3)) * t343 + (t561 * t373 + pkin(1) * t607 + ((-t137 * t367 + t161 * t373) * pkin(2) + (-0.2e1 * t137 * t564 + t161 * t699 - t161) * pkin(3)) * t354) * t355 - pkin(3) * t619 + t161 * t565 - t561 * t569 + pkin(3) * t137;
t628 = t62 * t165;
t534 = t322 * t628;
t636 = (t628 + (-t280 * t290 - pkin(3)) * t104) * t322;
t427 = (-t534 - t636) * t290;
t576 = t327 * t368;
t144 = t330 * t298 - t304 * t576;
t589 = t290 * t320;
t512 = t144 * t589;
t108 = t363 * t512;
t573 = t330 * t368;
t145 = t327 * t298 + t304 * t573;
t511 = t145 * t589;
t109 = t364 * t511;
t134 = t158 * t393 * t292 * t320;
t579 = t320 * t374;
t206 = t362 * t579;
t455 = t108 + t109 + t134 + t206;
t50 = t104 * t427 + t455;
t687 = t50 * pkin(1);
t369 = sin(qJ(3,1));
t568 = t369 * t354;
t166 = t375 * t355 - t568;
t163 = 0.1e1 / t166;
t560 = (-pkin(1) * t442 + t342 * t378) * t370 + (t342 * t442 + t338) * t376;
t563 = t375 * t369;
t596 = (-t369 * pkin(3) + t315) * t375;
t606 = t162 * t369;
t618 = t138 * t352;
t63 = ((t138 * t375 + t606) * pkin(2) + (0.2e1 * t162 * t563 - t138 + 0.2e1 * t618) * pkin(3)) * t343 + (t560 * t375 + pkin(1) * t606 + ((-t138 * t369 + t162 * t375) * pkin(2) + (-0.2e1 * t138 * t563 + t162 * t698 - t162) * pkin(3)) * t354) * t355 - pkin(3) * t618 + t162 * t596 - t560 * t568 + pkin(3) * t138;
t627 = t63 * t163;
t540 = t325 * t627;
t634 = (t627 + (-t280 * t293 - pkin(3)) * t105) * t325;
t429 = (-t540 - t634) * t293;
t575 = t328 * t370;
t146 = t331 * t299 - t305 * t575;
t586 = t293 * t323;
t510 = t146 * t586;
t110 = t363 * t510;
t572 = t331 * t370;
t147 = t328 * t299 + t305 * t572;
t509 = t147 * t586;
t111 = t364 * t509;
t135 = t159 * t393 * t295 * t323;
t578 = t323 * t376;
t207 = t362 * t578;
t454 = t110 + t111 + t135 + t207;
t51 = t105 * t429 + t454;
t686 = t51 * pkin(1);
t270 = t359 + t345;
t253 = qJ(1,3) + t270;
t685 = sin(t253) / 0.4e1;
t271 = -t359 + t345;
t254 = qJ(1,3) - t271;
t218 = sin(t254);
t684 = -t218 / 0.4e1;
t256 = qJ(1,3) - t270;
t220 = sin(t256);
t683 = -t220 / 0.4e1;
t272 = t360 + t346;
t257 = qJ(1,2) + t272;
t682 = sin(t257) / 0.4e1;
t273 = -t360 + t346;
t259 = qJ(1,2) - t273;
t223 = sin(t259);
t681 = -t223 / 0.4e1;
t260 = qJ(1,2) - t272;
t224 = sin(t260);
t680 = -t224 / 0.4e1;
t310 = qJ(1,1) + t347;
t261 = t361 + t310;
t679 = sin(t261) / 0.4e1;
t313 = qJ(1,1) - t347;
t263 = t361 + t313;
t227 = sin(t263);
t678 = -t227 / 0.4e1;
t264 = -t361 + t313;
t228 = sin(t264);
t677 = -t228 / 0.4e1;
t676 = -cos(t254) / 0.4e1;
t255 = qJ(1,3) + t271;
t675 = -cos(t255) / 0.4e1;
t258 = qJ(1,2) + t273;
t674 = -cos(t258) / 0.4e1;
t673 = -cos(t259) / 0.4e1;
t262 = -t361 + t310;
t672 = -cos(t262) / 0.4e1;
t671 = -cos(t263) / 0.4e1;
t269 = -t361 + t347;
t670 = sin(t269) / 0.2e1;
t669 = sin(t271) / 0.2e1;
t668 = sin(t273) / 0.2e1;
t667 = -0.2e1 * t343 + 0.1e1;
t666 = t701 - 0.2e1;
t311 = qJ(1,3) - t345;
t662 = g(3) * sin(t311);
t312 = qJ(1,2) - t346;
t661 = g(3) * sin(t312);
t660 = g(3) * sin(t313);
t659 = g(3) * cos(t311);
t658 = g(3) * cos(t312);
t657 = g(3) * cos(t313);
t653 = t350 * pkin(3);
t652 = t351 * pkin(3);
t651 = t352 * pkin(3);
t650 = t371 * pkin(2);
t649 = t373 * pkin(2);
t648 = t375 * pkin(2);
t477 = t370 * pkin(1) - t376 * t342;
t120 = t477 * t568 + (t352 - 0.1e1) * t370 * pkin(3);
t141 = pkin(1) * t369 + (-pkin(3) + t648 + 0.2e1 * t651) * t354;
t387 = -pkin(3) / 0.2e1;
t204 = t651 + t648 / 0.2e1 + t387;
t494 = t370 * t568;
t430 = pkin(2) * t494 + (t494 * t703 - t477) * t375;
t388 = pkin(2) / 0.2e1;
t593 = (t375 * pkin(3) + t388) * t369;
t80 = (-t204 * t575 + t331 * t593) * t702 + (t331 * t141 + t328 * t430) * t355 + t120 * t328 + t331 * t596;
t647 = t163 * t80;
t81 = (t204 * t572 + t328 * t593) * t702 + (t328 * t141 - t331 * t430) * t355 - t120 * t331 + t328 * t596;
t646 = t163 * t81;
t476 = pkin(1) * t366 - t372 * t340;
t118 = t476 * t570 + (t350 - 0.1e1) * t366 * pkin(3);
t139 = pkin(1) * t365 + (-pkin(3) + t650 + 0.2e1 * t653) * t354;
t202 = t653 + t650 / 0.2e1 + t387;
t496 = t366 * t570;
t432 = pkin(2) * t496 + (t496 * t703 - t476) * t371;
t595 = (t371 * pkin(3) + t388) * t365;
t76 = (-t202 * t577 + t329 * t595) * t702 + (t329 * t139 + t326 * t432) * t355 + t118 * t326 + t329 * t567;
t645 = t164 * t76;
t77 = (t202 * t574 + t326 * t595) * t702 + (t326 * t139 - t329 * t432) * t355 - t118 * t329 + t326 * t567;
t644 = t164 * t77;
t475 = pkin(1) * t368 - t374 * t341;
t119 = t475 * t569 + (t351 - 0.1e1) * t368 * pkin(3);
t140 = pkin(1) * t367 + (-pkin(3) + t649 + 0.2e1 * t652) * t354;
t203 = t652 + t649 / 0.2e1 + t387;
t495 = t368 * t569;
t431 = pkin(2) * t495 + (t495 * t703 - t475) * t373;
t594 = (t373 * pkin(3) + t388) * t367;
t78 = (-t203 * t576 + t330 * t594) * t702 + (t330 * t140 + t327 * t431) * t355 + t119 * t327 + t330 * t565;
t643 = t165 * t78;
t79 = (t203 * t573 + t327 * t594) * t702 + (t327 * t140 - t330 * t431) * t355 - t119 * t330 + t327 * t565;
t642 = t165 * t79;
t97 = t103 ^ 2;
t641 = t288 * t97;
t98 = t104 ^ 2;
t640 = t291 * t98;
t99 = t105 ^ 2;
t639 = t294 * t99;
t637 = t318 * t97;
t635 = t321 * t98;
t633 = t324 * t99;
t632 = t350 * t49;
t631 = t351 * t50;
t630 = t352 * t51;
t626 = t103 * t287;
t625 = t104 * t290;
t624 = t105 * t293;
t130 = t172 * t372 + t366 * t340;
t623 = t130 * t317;
t131 = t173 * t374 + t368 * t341;
t622 = t131 * t320;
t132 = t174 * t376 + t370 * t342;
t621 = t132 * t323;
t617 = t142 * t287;
t616 = t143 * t287;
t615 = t144 * t290;
t614 = t145 * t290;
t613 = t146 * t293;
t612 = t147 * t293;
t605 = t163 * t323;
t604 = t164 * t317;
t603 = t165 * t320;
t601 = t175 * t372;
t599 = t176 * t374;
t597 = t177 * t376;
t591 = t287 * t326;
t590 = t287 * t329;
t588 = t290 * t327;
t587 = t290 * t330;
t585 = t293 * t328;
t584 = t293 * t331;
t571 = t354 * t355;
t386 = 0.2e1 * pkin(7);
t344 = t386 + qJ(3,2);
t296 = sin(t344);
t559 = -t296 - t367;
t348 = qJ(3,3) + t386;
t300 = sin(t348);
t558 = -t300 - t365;
t349 = qJ(3,1) + t386;
t301 = sin(t349);
t557 = -t301 - t369;
t302 = cos(t344);
t556 = t302 + t373;
t306 = cos(t348);
t555 = t306 + t371;
t307 = cos(t349);
t554 = t307 + t375;
t552 = 0.2e1 * t688;
t551 = 0.2e1 * t687;
t550 = 0.2e1 * t686;
t549 = t343 - 0.1e1 / 0.2e1;
t548 = t350 - 0.1e1 / 0.2e1;
t547 = t351 - 0.1e1 / 0.2e1;
t546 = t352 - 0.1e1 / 0.2e1;
t545 = 2 * t393;
t544 = -0.8e1 * t571;
t543 = 0.4e1 * t571;
t542 = t80 * t605;
t541 = t81 * t605;
t539 = t76 * t604;
t538 = t77 * t604;
t536 = t78 * t603;
t535 = t79 * t603;
t533 = t49 * t591;
t532 = t49 * t590;
t531 = t288 * t637;
t530 = t319 * t641;
t529 = t289 * t637;
t528 = t50 * t588;
t527 = t50 * t587;
t526 = t291 * t635;
t525 = t322 * t640;
t524 = t292 * t635;
t523 = t51 * t585;
t522 = t51 * t584;
t521 = t294 * t633;
t520 = t325 * t639;
t519 = t295 * t633;
t518 = t49 * t566;
t517 = t50 * t564;
t516 = t51 * t563;
t515 = pkin(2) * t703;
t333 = g(3) * t366;
t151 = t333 - t601;
t508 = t151 * t592;
t334 = g(3) * t368;
t152 = t334 - t599;
t507 = t152 * t589;
t335 = g(3) * t370;
t153 = t335 - t597;
t506 = t153 * t586;
t505 = t154 * t592;
t504 = t155 * t589;
t503 = t156 * t586;
t502 = t288 * t611;
t501 = t291 * t610;
t500 = t294 * t609;
t492 = t684 - cos(t271) / 0.2e1;
t491 = t220 / 0.4e1 - cos(t270) / 0.2e1;
t490 = t681 - cos(t273) / 0.2e1;
t489 = t224 / 0.4e1 - cos(t272) / 0.2e1;
t488 = t678 - cos(t269) / 0.2e1;
t268 = t361 + t347;
t487 = t228 / 0.4e1 - cos(t268) / 0.2e1;
t232 = cos(t256);
t486 = -t232 / 0.4e1 + t676;
t485 = t232 / 0.4e1 + t676;
t236 = cos(t260);
t484 = -t236 / 0.4e1 + t673;
t483 = t236 / 0.4e1 + t673;
t240 = cos(t264);
t482 = -t240 / 0.4e1 + t671;
t481 = t240 / 0.4e1 + t671;
t480 = pkin(2) * t531;
t479 = pkin(2) * t526;
t478 = pkin(2) * t521;
t474 = t626 * t629;
t473 = t625 * t628;
t472 = t624 * t627;
t471 = t163 * t520;
t470 = t164 * t530;
t469 = t165 * t525;
t468 = t326 * t529;
t467 = t329 * t529;
t466 = t327 * t524;
t465 = t330 * t524;
t464 = t328 * t519;
t463 = t331 * t519;
t459 = t566 * t571;
t458 = t564 * t571;
t457 = t563 * t571;
t219 = sin(t255);
t229 = cos(t253);
t308 = qJ(1,3) + t345;
t450 = g(1) * t675 + g(2) * t685 + t219 * t694 + t229 * t696 + sin(t308) * t692;
t222 = sin(t258);
t233 = cos(t257);
t309 = qJ(1,2) + t346;
t449 = g(1) * t674 + g(2) * t682 + t222 * t694 + t233 * t696 + sin(t309) * t692;
t226 = sin(t262);
t237 = cos(t261);
t448 = g(1) * t672 + g(2) * t679 + t226 * t694 + t237 * t696 + sin(t310) * t692;
t447 = g(1) * t685 + g(2) * t675 + t219 * t695 + t229 * t693 + cos(t308) * t692;
t446 = g(1) * t682 + g(2) * t674 + t222 * t695 + t233 * t693 + cos(t309) * t692;
t445 = g(1) * t679 + g(2) * t672 + t226 * t695 + t237 * t693 + cos(t310) * t692;
t281 = 0.2e1 * t345;
t389 = (qJ(2,3) ^ 2);
t392 = pkin(3) ^ 2;
t395 = pkin(2) ^ 2;
t423 = -t395 * cos(t386) - (2 * pkin(1) ^ 2) - (2 * pkin(6) ^ 2) - t392 - t395 + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t441 = t362 * t623 + t363 * t539 + t364 * t538 + (-t103 * sin(t281) / 0.2e1 + t160 * t393 * t172) * t287 * t499 - (t61 * t604 * t704 + 0.4e1 * (-t314 - t665) * (pkin(1) * t626 - t629 / 0.2e1) * t317 + (0.2e1 * t340 * t583 + ((qJ(2,3) * t697) - t392 * cos(t281) - (2 * t389) - t555 * t515 + t423) * t317 * t103) * t287) * t318 * t626 / 0.2e1 - t172 * t319 * t474;
t282 = 0.2e1 * t346;
t390 = (qJ(2,2) ^ 2);
t440 = t362 * t622 + t363 * t536 + t364 * t535 + (-t104 * sin(t282) / 0.2e1 + t161 * t393 * t173) * t290 * t498 - (t62 * t603 * t704 + 0.4e1 * (-t314 - t664) * (pkin(1) * t625 - t628 / 0.2e1) * t320 + (0.2e1 * t341 * t582 + ((qJ(2,2) * t697) - t392 * cos(t282) - (2 * t390) - t556 * t515 + t423) * t320 * t104) * t290) * t321 * t625 / 0.2e1 - t173 * t322 * t473;
t283 = 0.2e1 * t347;
t391 = (qJ(2,1) ^ 2);
t439 = t362 * t621 + t363 * t542 + t364 * t541 + (-t105 * sin(t283) / 0.2e1 + t162 * t393 * t174) * t293 * t497 - (t63 * t605 * t704 + 0.4e1 * (-t314 - t663) * (pkin(1) * t624 - t627 / 0.2e1) * t323 + (0.2e1 * t342 * t581 + ((qJ(2,1) * t697) - t392 * cos(t283) - (2 * t391) - t554 * t515 + t423) * t323 * t105) * t293) * t324 * t624 / 0.2e1 - t174 * t325 * t472;
t438 = 0.2e1 * t453;
t437 = 0.2e1 * t452;
t436 = 0.2e1 * t451;
t435 = t462 * t566;
t434 = t461 * t564;
t433 = t460 * t563;
t426 = -t441 + t552;
t425 = -t440 + t551;
t424 = -t439 + t550;
t422 = t354 * t49 + t355 * t438;
t421 = t354 * t438 - t355 * t49;
t420 = t354 * t50 + t355 * t437;
t419 = t354 * t437 - t355 * t50;
t418 = t354 * t51 + t355 * t436;
t417 = t354 * t436 - t355 * t51;
t416 = -t441 - t601;
t415 = -t440 - t599;
t414 = -t439 - t597;
t413 = t356 * t502 - t426;
t412 = t357 * t501 - t425;
t411 = t358 * t500 - t424;
t410 = 0.2e1 * t318 * t474 - t154;
t409 = 0.2e1 * t321 * t473 - t155;
t408 = 0.2e1 * t324 * t472 - t156;
t382 = pkin(1) * g(3);
t245 = sin(t272);
t243 = sin(t270);
t241 = sin(t268);
t171 = t375 * t354 + t369 * t355;
t170 = t373 * t354 + t367 * t355;
t169 = t371 * t354 + t365 * t355;
t129 = t666 * t352 - 0.4e1 * t457 + t667;
t128 = t666 * t351 - 0.4e1 * t458 + t667;
t127 = t666 * t350 - 0.4e1 * t459 + t667;
t123 = t546 * t571 + t549 * t563;
t122 = t547 * t571 + t549 * t564;
t121 = t548 * t571 + t549 * t566;
t96 = t102 * t354 + t355 * t500;
t95 = t101 * t354 + t355 * t501;
t94 = t100 * t354 + t355 * t502;
t93 = -t102 * t355 + t354 * t500;
t92 = -t101 * t355 + t354 * t501;
t91 = -t100 * t355 + t354 * t502;
t75 = -t96 * t369 - t93 * t375;
t74 = -t95 * t367 - t92 * t373;
t73 = -t94 * t365 - t91 * t371;
t72 = -t93 * t369 + t96 * t375;
t71 = -t92 * t367 + t95 * t373;
t70 = -t91 * t365 + t94 * t371;
t48 = 0.2e1 * qJ(2,1) * t51 + t408;
t47 = 0.2e1 * qJ(2,2) * t50 + t409;
t46 = 0.2e1 * qJ(2,3) * t49 + t410;
t42 = -t51 * t358 + (pkin(1) * t639 - 0.2e1 * t472) * t324;
t41 = -t50 * t357 + (pkin(1) * t640 - 0.2e1 * t473) * t321;
t40 = -t49 * t356 + (pkin(1) * t641 - 0.2e1 * t474) * t318;
t36 = -t660 / 0.2e1 + t42 * t299 + (t301 / 0.2e1 + t369 / 0.2e1) * t478 + (t487 + t488) * g(2) + (t670 - t241 / 0.2e1 - t482) * g(1) + t448;
t35 = -t661 / 0.2e1 + t41 * t298 + (t296 / 0.2e1 + t367 / 0.2e1) * t479 + (t489 + t490) * g(2) + (t668 - t245 / 0.2e1 - t484) * g(1) + t449;
t34 = -t662 / 0.2e1 + t40 * t297 + (t300 / 0.2e1 + t365 / 0.2e1) * t480 + (t491 + t492) * g(2) + (t669 - t243 / 0.2e1 - t486) * g(1) + t450;
t33 = t657 / 0.2e1 + t42 * t305 + (t307 / 0.2e1 + t375 / 0.2e1) * t478 + (t670 + t241 / 0.2e1 - t481) * g(2) + (t487 - t488) * g(1) + t445;
t32 = t658 / 0.2e1 + t41 * t304 + (t302 / 0.2e1 + t373 / 0.2e1) * t479 + (t668 + t245 / 0.2e1 - t483) * g(2) + (t489 - t490) * g(1) + t446;
t31 = t659 / 0.2e1 + t40 * t303 + (t306 / 0.2e1 + t371 / 0.2e1) * t480 + (t669 + t243 / 0.2e1 - t485) * g(2) + (t491 - t492) * g(1) + t447;
t30 = -t369 * t417 + t375 * t418;
t29 = -t367 * t419 + t373 * t420;
t28 = -t365 * t421 + t371 * t422;
t27 = t369 * t418 + t375 * t417;
t26 = t367 * t420 + t373 * t419;
t25 = t365 * t422 + t371 * t421;
t24 = t335 + t414 + t550;
t23 = t334 + t415 + t551;
t22 = t333 + t416 + t552;
t21 = (-t153 - t424) * t354;
t20 = (-t152 - t425) * t354;
t19 = (-t151 - t426) * t354;
t18 = -qJ(2,1) * t521 - t153 + t439 - t686;
t17 = -qJ(2,2) * t526 - t152 + t440 - t687;
t16 = -qJ(2,3) * t531 - t151 + t441 - t688;
t15 = (-0.2e1 * t630 + (t429 + 0.4e1 * t433) * t105 + t454) * t343 + (t516 / 0.2e1 + t546 * t451) * t543 - 0.2e1 * t105 * t433 + t630;
t14 = (-0.2e1 * t631 + (t427 + 0.4e1 * t434) * t104 + t455) * t343 + (t517 / 0.2e1 + t547 * t452) * t543 - 0.2e1 * t104 * t434 + t631;
t13 = (-0.2e1 * t632 + (t428 + 0.4e1 * t435) * t103 + t456) * t343 + (t518 / 0.2e1 + t548 * t453) * t543 - 0.2e1 * t103 * t435 + t632;
t12 = (-t630 / 0.2e1 + t111 / 0.4e1 + t110 / 0.4e1 + t207 / 0.4e1 + t135 / 0.4e1) * t544 + t666 * t516 + ((-t634 / 0.4e1 - t540 / 0.4e1) * t293 * t544 + ((t698 - 0.1e1) * t701 - 0.8e1 * t457 - 0.4e1 * t352 + 0.2e1) * t460) * t105;
t11 = (-t631 / 0.2e1 + t109 / 0.4e1 + t108 / 0.4e1 + t206 / 0.4e1 + t134 / 0.4e1) * t544 + t666 * t517 + ((-t636 / 0.4e1 - t534 / 0.4e1) * t290 * t544 + ((t699 - 0.1e1) * t701 - 0.8e1 * t458 - 0.4e1 * t351 + 0.2e1) * t461) * t104;
t10 = (-t632 / 0.2e1 + t107 / 0.4e1 + t106 / 0.4e1 + t205 / 0.4e1 + t133 / 0.4e1) * t544 + t666 * t518 + ((-t638 / 0.4e1 - t537 / 0.4e1) * t287 * t544 + ((t700 - 0.1e1) * t701 - 0.8e1 * t459 - 0.4e1 * t350 + 0.2e1) * t462) * t103;
t9 = t382 * t370 + t51 * t391 + t408 * qJ(2,1) + (t414 + t686) * pkin(1);
t8 = t382 * t368 + t50 * t390 + t409 * qJ(2,2) + (t415 + t687) * pkin(1);
t7 = t382 * t366 + t49 * t389 + t410 * qJ(2,3) + (t416 + t688) * pkin(1);
t6 = t660 / 0.2e1 - t411 * t305 + t299 * t689 + (t677 + t227 / 0.4e1) * g(2) + t482 * g(1) + (t451 * t557 + t51 * t554) * pkin(2) + t448;
t5 = t661 / 0.2e1 - t412 * t304 + t298 * t690 + (t680 + t223 / 0.4e1) * g(2) + t484 * g(1) + (t452 * t559 + t556 * t50) * pkin(2) + t449;
t4 = t662 / 0.2e1 - t413 * t303 + t297 * t691 + (t683 + t218 / 0.4e1) * g(2) + t486 * g(1) + (t453 * t558 + t49 * t555) * pkin(2) + t450;
t3 = -t657 / 0.2e1 + t411 * t299 + t305 * t689 + t481 * g(2) + (t677 + t678) * g(1) + (-t451 * t554 + t51 * t557) * pkin(2) + t445;
t2 = -t658 / 0.2e1 + t412 * t298 + t304 * t690 + t483 * g(2) + (t680 + t681) * g(1) + (-t556 * t452 + t50 * t559) * pkin(2) + t446;
t1 = -t659 / 0.2e1 + t413 * t297 + t303 * t691 + t485 * g(2) + (t683 + t684) * g(1) + (-t453 * t555 + t49 * t558) * pkin(2) + t447;
t37 = [t49 * t513 + t50 * t511 + t509 * t51, t143 * t508 + t145 * t507 + t147 * t506, t143 * t505 + t145 * t504 + t147 * t503, ((t24 * t612 - t51 * t646) * t323 + (t23 * t614 - t50 * t642) * t320 + (t22 * t616 - t49 * t644) * t317) * t355, t19 * t513 + t20 * t511 + t21 * t509 + (t49 * t538 + t50 * t535 + t51 * t541) * t354, t46 * t513 - t469 * t79 + t47 * t511 - t470 * t77 - t471 * t81 + t48 * t509, (t18 * t646 + t612 * t9) * t323 + (t17 * t642 + t614 * t8) * t320 + (t16 * t644 + t616 * t7) * t317, t13 * t513 + t14 * t511 + t15 * t509 + (-t121 * t468 - t122 * t466 - t123 * t464) * t545, t10 * t513 + t11 * t511 + t12 * t509 + (-t127 * t468 - t128 * t466 - t129 * t464) * t393, t70 * t513 + t71 * t511 + t72 * t509 + (t169 * t533 + t170 * t528 + t171 * t523) * t393, t73 * t513 + t74 * t511 + t75 * t509 + (t166 * t523 + t167 * t533 + t168 * t528) * t393, (t100 * t591 + t101 * t588 + t102 * t585) * t393, (t27 * t646 + t6 * t612) * t323 + (t26 * t642 + t5 * t614) * t320 + (t25 * t644 + t4 * t616) * t317 + (t34 * t591 + t35 * t588 + t36 * t585) * t393, (t3 * t612 + t30 * t646) * t323 + (t2 * t614 + t29 * t642) * t320 + (t1 * t616 + t28 * t644) * t317 + (t31 * t591 + t32 * t588 + t33 * t585) * t393, t364 - g(1); t49 * t514 + t50 * t512 + t51 * t510, t142 * t508 + t144 * t507 + t146 * t506, t142 * t505 + t144 * t504 + t146 * t503, ((t24 * t613 - t51 * t647) * t323 + (t23 * t615 - t50 * t643) * t320 + (t22 * t617 - t49 * t645) * t317) * t355, t19 * t514 + t20 * t512 + t21 * t510 + (t49 * t539 + t50 * t536 + t51 * t542) * t354, t46 * t514 - t469 * t78 + t47 * t512 - t470 * t76 - t471 * t80 + t48 * t510, (t18 * t647 + t613 * t9) * t323 + (t17 * t643 + t615 * t8) * t320 + (t16 * t645 + t617 * t7) * t317, t13 * t514 + t14 * t512 + t15 * t510 + (-t121 * t467 - t122 * t465 - t123 * t463) * t545, t10 * t514 + t11 * t512 + t12 * t510 + (-t127 * t467 - t128 * t465 - t129 * t463) * t393, t70 * t514 + t71 * t512 + t72 * t510 + (t169 * t532 + t170 * t527 + t171 * t522) * t393, t73 * t514 + t74 * t512 + t75 * t510 + (t166 * t522 + t167 * t532 + t168 * t527) * t393, (t100 * t590 + t101 * t587 + t102 * t584) * t393, (t27 * t647 + t6 * t613) * t323 + (t26 * t643 + t5 * t615) * t320 + (t25 * t645 + t4 * t617) * t317 + (t34 * t590 + t35 * t587 + t36 * t584) * t393, (t3 * t613 + t30 * t647) * t323 + (t2 * t615 + t29 * t643) * t320 + (t1 * t617 + t28 * t645) * t317 + (t31 * t590 + t32 * t587 + t33 * t584) * t393, t363 - g(2); t49 * t580 + t50 * t579 + t51 * t578, t151 * t580 + t152 * t579 + t153 * t578, t154 * t580 + t155 * t579 + t156 * t578, ((-t132 * t51 + t24 * t376) * t323 + (-t131 * t50 + t23 * t374) * t320 + (-t130 * t49 + t22 * t372) * t317) * t355, t19 * t580 + t20 * t579 + t21 * t578 + (t49 * t623 + t50 * t622 + t51 * t621) * t354, -t130 * t530 - t131 * t525 - t132 * t520 + t46 * t580 + t47 * t579 + t48 * t578, (t132 * t18 + t376 * t9) * t323 + (t131 * t17 + t374 * t8) * t320 + (t130 * t16 + t372 * t7) * t317, t13 * t580 + t14 * t579 + t15 * t578, t10 * t580 + t11 * t579 + t12 * t578, t578 * t72 + t579 * t71 + t580 * t70, t578 * t75 + t579 * t74 + t580 * t73, 0, (t132 * t27 + t376 * t6) * t323 + (t131 * t26 + t374 * t5) * t320 + (t130 * t25 + t372 * t4) * t317, (t132 * t30 + t3 * t376) * t323 + (t131 * t29 + t2 * t374) * t320 + (t1 * t372 + t130 * t28) * t317, t362 - g(3);];
tauX_reg  = t37;
