% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR12V2G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRPRR12V2G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:21:12
% EndTime: 2020-08-06 19:21:43
% DurationCPUTime: 31.09s
% Computational Cost: add. (263400->887), mult. (352557->1553), div. (15066->12), fcn. (187377->18), ass. (0->558)
t692 = 2 * pkin(1);
t315 = sin(qJ(2,3));
t270 = pkin(1) * t315 + qJ(3,3);
t309 = legFrame(3,2);
t282 = sin(t309);
t285 = cos(t309);
t321 = cos(qJ(2,3));
t304 = t321 ^ 2;
t316 = sin(qJ(1,3));
t322 = cos(qJ(1,3));
t327 = xDP(3);
t328 = xDP(2);
t329 = xDP(1);
t381 = t282 * t329 + t285 * t328;
t279 = t315 * qJ(3,3);
t267 = t279 + pkin(1);
t330 = pkin(5) - pkin(6);
t520 = t330 * t322;
t384 = t267 * t316 - t520;
t276 = t316 * t330;
t498 = t267 * t322 + t276;
t331 = pkin(2) + pkin(3);
t513 = t331 * t315;
t521 = t329 * t331;
t522 = t328 * t331;
t523 = t327 * t331;
t160 = ((-qJ(3,3) * t328 + t316 * t521) * t285 + (-qJ(3,3) * t329 - t316 * t522) * t282 + t322 * t523) * t304 + ((t328 * t513 + t329 * t384) * t285 + (-t328 * t384 + t329 * t513) * t282 + t498 * t327) * t321 + t381 * t270;
t380 = -t282 * t328 + t285 * t329;
t202 = -t316 * t327 + t322 * t380;
t510 = t331 * t321;
t244 = t510 + t267;
t230 = 0.1e1 / t244 ^ 2;
t334 = 0.1e1 / qJ(3,3);
t587 = t230 * t334;
t452 = t202 * t587;
t406 = t160 * t452;
t229 = 0.1e1 / t244;
t333 = qJ(3,3) ^ 2;
t297 = pkin(1) ^ 2 + pkin(5) ^ 2;
t400 = t297 + (-2 * pkin(5) + pkin(6)) * pkin(6);
t250 = t333 + t400;
t490 = 0.2e1 * t331;
t265 = pkin(1) * t522;
t266 = pkin(1) * t521;
t415 = pkin(1) * t316 - t520;
t357 = (pkin(1) * t322 + t276) * t327 + t380 * t415;
t547 = (qJ(3,3) + t331) * (-qJ(3,3) + t331);
t428 = t322 * t547;
t429 = t316 * t547;
t485 = -0.2e1 * t521;
t486 = -0.2e1 * t522;
t526 = t327 * t322;
t542 = t316 * qJ(3,3);
t430 = t315 * t547;
t680 = -pkin(1) * qJ(3,3) + t430;
t701 = 0.2e1 * t279;
t139 = ((qJ(3,3) * t486 + t329 * t429) * t285 + (qJ(3,3) * t485 - t328 * t429) * t282 + t327 * t428) * t304 + (((t316 * t380 + t526) * t701 + t357) * t331 + t680 * t381) * t321 + (qJ(3,3) * t357 + t265 * t285 + t266 * t282) * t315 - qJ(3,3) * ((-t329 * t542 - t522) * t285 + (t328 * t542 - t521) * t282 - qJ(3,3) * t526);
t588 = t229 * t334;
t136 = t139 * t588;
t463 = t160 * t588;
t154 = pkin(2) * t463;
t122 = t154 - t136;
t409 = pkin(3) * t463;
t115 = -t409 - t122;
t112 = t330 * t115;
t607 = t202 * t229;
t479 = qJ(3,3) * t607;
t504 = t479 * t692 + t112;
t622 = t160 * t330;
t103 = t504 * t315 + (t321 * t622 + (t267 * t321 * t490 + t304 * t547 + t250) * t202) * t229;
t420 = t315 * t542;
t218 = t415 + t420;
t541 = t316 * t331;
t559 = t285 * qJ(3,3);
t181 = (-t282 * t541 - t559) * t304 + (-t218 * t282 + t285 * t513) * t321 + t285 * t270;
t571 = t282 * qJ(3,3);
t182 = (t285 * t541 - t571) * t304 + (t218 * t285 + t282 * t513) * t321 + t282 * t270;
t312 = xDDP(3);
t313 = xDDP(2);
t314 = xDDP(1);
t205 = t322 * t510 + t498;
t595 = t205 * t321;
t360 = t181 * t313 + t182 * t314 + t312 * t595;
t605 = t202 * t330;
t451 = t315 * t605;
t507 = t331 * t334;
t145 = (-t160 * t507 + t451) * t229;
t454 = t229 * t605;
t626 = qJ(3,3) * t304;
t387 = (-t145 * t321 + t267 * t463) * t139 + t160 * (-t454 * t626 + ((-pkin(3) * t160 * t334 + t451) * t229 - t122) * t510 + t115 * t267);
t335 = 0.1e1 / qJ(3,3) ^ 2;
t586 = t230 * t335;
t606 = t202 * t230;
t46 = -t387 * t586 + (t103 * t321 * t606 + t229 * t360) * t334;
t34 = pkin(5) * t46 + t406 * t692;
t707 = t34 * t315;
t317 = sin(qJ(2,2));
t271 = pkin(1) * t317 + qJ(3,2);
t310 = legFrame(2,2);
t283 = sin(t310);
t286 = cos(t310);
t323 = cos(qJ(2,2));
t305 = t323 ^ 2;
t318 = sin(qJ(1,2));
t324 = cos(qJ(1,2));
t379 = t283 * t329 + t286 * t328;
t280 = t317 * qJ(3,2);
t268 = t280 + pkin(1);
t519 = t330 * t324;
t383 = t268 * t318 - t519;
t277 = t318 * t330;
t497 = t268 * t324 + t277;
t512 = t331 * t317;
t161 = ((-qJ(3,2) * t328 + t318 * t521) * t286 + (-qJ(3,2) * t329 - t318 * t522) * t283 + t324 * t523) * t305 + ((t328 * t512 + t329 * t383) * t286 + (-t328 * t383 + t329 * t512) * t283 + t497 * t327) * t323 + t379 * t271;
t378 = -t283 * t328 + t286 * t329;
t203 = -t318 * t327 + t324 * t378;
t509 = t331 * t323;
t245 = t509 + t268;
t233 = 0.1e1 / t245 ^ 2;
t337 = 0.1e1 / qJ(3,2);
t583 = t233 * t337;
t448 = t203 * t583;
t405 = t161 * t448;
t232 = 0.1e1 / t245;
t336 = qJ(3,2) ^ 2;
t251 = t336 + t400;
t414 = pkin(1) * t318 - t519;
t356 = (pkin(1) * t324 + t277) * t327 + t378 * t414;
t546 = (qJ(3,2) + t331) * (-qJ(3,2) + t331);
t425 = t324 * t546;
t426 = t318 * t546;
t525 = t327 * t324;
t538 = t318 * qJ(3,2);
t427 = t317 * t546;
t681 = -pkin(1) * qJ(3,2) + t427;
t700 = 0.2e1 * t280;
t140 = ((qJ(3,2) * t486 + t329 * t426) * t286 + (qJ(3,2) * t485 - t328 * t426) * t283 + t327 * t425) * t305 + (((t318 * t378 + t525) * t700 + t356) * t331 + t681 * t379) * t323 + (qJ(3,2) * t356 + t265 * t286 + t266 * t283) * t317 - qJ(3,2) * ((-t329 * t538 - t522) * t286 + (t328 * t538 - t521) * t283 - qJ(3,2) * t525);
t584 = t232 * t337;
t137 = t140 * t584;
t461 = t161 * t584;
t155 = pkin(2) * t461;
t124 = t155 - t137;
t408 = pkin(3) * t461;
t116 = -t408 - t124;
t113 = t330 * t116;
t603 = t203 * t232;
t480 = qJ(3,2) * t603;
t503 = t480 * t692 + t113;
t621 = t161 * t330;
t104 = t503 * t317 + (t323 * t621 + (t268 * t323 * t490 + t305 * t546 + t251) * t203) * t232;
t418 = t317 * t538;
t220 = t414 + t418;
t537 = t318 * t331;
t555 = t286 * qJ(3,2);
t183 = (-t283 * t537 - t555) * t305 + (-t220 * t283 + t286 * t512) * t323 + t286 * t271;
t567 = t283 * qJ(3,2);
t184 = (t286 * t537 - t567) * t305 + (t220 * t286 + t283 * t512) * t323 + t283 * t271;
t206 = t324 * t509 + t497;
t593 = t206 * t323;
t359 = t183 * t313 + t184 * t314 + t312 * t593;
t601 = t203 * t330;
t447 = t317 * t601;
t506 = t331 * t337;
t146 = (-t161 * t506 + t447) * t232;
t450 = t232 * t601;
t628 = qJ(3,2) * t305;
t386 = (-t146 * t323 + t268 * t461) * t140 + t161 * (-t450 * t628 + ((-pkin(3) * t161 * t337 + t447) * t232 - t124) * t509 + t116 * t268);
t338 = 0.1e1 / qJ(3,2) ^ 2;
t582 = t233 * t338;
t602 = t203 * t233;
t47 = -t386 * t582 + (t104 * t323 * t602 + t232 * t359) * t337;
t35 = t47 * pkin(5) + t405 * t692;
t706 = t35 * t317;
t319 = sin(qJ(2,1));
t272 = pkin(1) * t319 + qJ(3,1);
t311 = legFrame(1,2);
t284 = sin(t311);
t287 = cos(t311);
t325 = cos(qJ(2,1));
t306 = t325 ^ 2;
t320 = sin(qJ(1,1));
t326 = cos(qJ(1,1));
t377 = t284 * t329 + t287 * t328;
t281 = t319 * qJ(3,1);
t269 = t281 + pkin(1);
t518 = t330 * t326;
t382 = t269 * t320 - t518;
t278 = t320 * t330;
t496 = t269 * t326 + t278;
t511 = t331 * t319;
t162 = ((-qJ(3,1) * t328 + t320 * t521) * t287 + (-qJ(3,1) * t329 - t320 * t522) * t284 + t326 * t523) * t306 + ((t328 * t511 + t329 * t382) * t287 + (-t328 * t382 + t329 * t511) * t284 + t496 * t327) * t325 + t377 * t272;
t376 = -t284 * t328 + t287 * t329;
t204 = -t320 * t327 + t326 * t376;
t508 = t331 * t325;
t246 = t508 + t269;
t236 = 0.1e1 / t246 ^ 2;
t340 = 0.1e1 / qJ(3,1);
t579 = t236 * t340;
t444 = t204 * t579;
t404 = t162 * t444;
t235 = 0.1e1 / t246;
t339 = qJ(3,1) ^ 2;
t252 = t339 + t400;
t413 = pkin(1) * t320 - t518;
t355 = (pkin(1) * t326 + t278) * t327 + t376 * t413;
t545 = (qJ(3,1) + t331) * (-qJ(3,1) + t331);
t422 = t326 * t545;
t423 = t320 * t545;
t524 = t327 * t326;
t534 = t320 * qJ(3,1);
t424 = t319 * t545;
t682 = -pkin(1) * qJ(3,1) + t424;
t699 = 0.2e1 * t281;
t141 = ((qJ(3,1) * t486 + t329 * t423) * t287 + (qJ(3,1) * t485 - t328 * t423) * t284 + t327 * t422) * t306 + (((t320 * t376 + t524) * t699 + t355) * t331 + t682 * t377) * t325 + (qJ(3,1) * t355 + t265 * t287 + t266 * t284) * t319 - qJ(3,1) * ((-t329 * t534 - t522) * t287 + (t328 * t534 - t521) * t284 - qJ(3,1) * t524);
t580 = t235 * t340;
t138 = t141 * t580;
t459 = t162 * t580;
t156 = pkin(2) * t459;
t126 = t156 - t138;
t407 = pkin(3) * t459;
t117 = -t407 - t126;
t114 = t330 * t117;
t599 = t204 * t235;
t481 = qJ(3,1) * t599;
t502 = t481 * t692 + t114;
t620 = t162 * t330;
t105 = t502 * t319 + (t325 * t620 + (t269 * t325 * t490 + t306 * t545 + t252) * t204) * t235;
t416 = t319 * t534;
t222 = t413 + t416;
t533 = t320 * t331;
t551 = t287 * qJ(3,1);
t185 = (-t284 * t533 - t551) * t306 + (-t222 * t284 + t287 * t511) * t325 + t287 * t272;
t563 = t284 * qJ(3,1);
t186 = (t287 * t533 - t563) * t306 + (t222 * t287 + t284 * t511) * t325 + t284 * t272;
t207 = t326 * t508 + t496;
t591 = t207 * t325;
t358 = t185 * t313 + t186 * t314 + t312 * t591;
t597 = t204 * t330;
t443 = t319 * t597;
t505 = t331 * t340;
t147 = (-t162 * t505 + t443) * t235;
t446 = t235 * t597;
t630 = qJ(3,1) * t306;
t385 = (-t147 * t325 + t269 * t459) * t141 + t162 * (-t446 * t630 + ((-pkin(3) * t162 * t340 + t443) * t235 - t126) * t508 + t117 * t269);
t341 = 0.1e1 / qJ(3,1) ^ 2;
t578 = t236 * t341;
t598 = t204 * t236;
t48 = -t385 * t578 + (t105 * t325 * t598 + t235 * t358) * t340;
t36 = t48 * pkin(5) + t404 * t692;
t705 = t36 * t319;
t237 = t235 * t236;
t535 = t319 * t340;
t431 = t237 * t535;
t550 = t287 * t326;
t432 = t235 * t550;
t562 = t284 * t326;
t433 = t235 * t562;
t528 = t325 * qJ(3,1);
t581 = t235 * t320;
t596 = t204 * t340;
t102 = t314 * t432 - t313 * t433 - t312 * t581 - (t138 * t319 + (t597 + (-t319 * t505 + t325) * t162) * t235) * t598 + t237 * (t511 - t528) * t162 * t596 - t204 * t141 * t431;
t243 = t287 * g(1) - t284 * g(2);
t210 = g(3) * t320 - t243 * t326;
t159 = t162 ^ 2;
t464 = t159 * t578;
t696 = -pkin(5) * t464 + t102 * t692 + t210;
t704 = t696 * t319;
t234 = t232 * t233;
t539 = t317 * t337;
t434 = t234 * t539;
t554 = t286 * t324;
t435 = t232 * t554;
t566 = t283 * t324;
t436 = t232 * t566;
t530 = t323 * qJ(3,2);
t585 = t232 * t318;
t600 = t203 * t337;
t101 = t314 * t435 - t313 * t436 - t312 * t585 - (t137 * t317 + (t601 + (-t317 * t506 + t323) * t161) * t232) * t602 + t234 * (t512 - t530) * t161 * t600 - t203 * t140 * t434;
t242 = t286 * g(1) - t283 * g(2);
t209 = g(3) * t318 - t242 * t324;
t158 = t161 ^ 2;
t465 = t158 * t582;
t697 = -pkin(5) * t465 + t101 * t692 + t209;
t703 = t697 * t317;
t231 = t229 * t230;
t543 = t315 * t334;
t437 = t231 * t543;
t558 = t285 * t322;
t438 = t229 * t558;
t570 = t282 * t322;
t439 = t229 * t570;
t532 = t321 * qJ(3,3);
t589 = t229 * t316;
t604 = t202 * t334;
t100 = t314 * t438 - t313 * t439 - t312 * t589 - (t136 * t315 + (t605 + (-t315 * t507 + t321) * t160) * t229) * t606 + t231 * (t513 - t532) * t160 * t604 - t202 * t139 * t437;
t241 = t285 * g(1) - t282 * g(2);
t208 = g(3) * t316 - t241 * t322;
t157 = t160 ^ 2;
t466 = t157 * t586;
t698 = -pkin(5) * t466 + t100 * t692 + t208;
t702 = t698 * t315;
t695 = t241 * t316;
t694 = t242 * t318;
t693 = t243 * t320;
t691 = 0.4e1 * t304;
t690 = 0.4e1 * t305;
t689 = 0.4e1 * t306;
t688 = 0.2e1 * t315;
t687 = 0.2e1 * t317;
t686 = 0.2e1 * t319;
t678 = 0.2e1 * pkin(2);
t674 = -0.2e1 * qJ(3,1);
t673 = -0.2e1 * qJ(3,2);
t672 = -0.2e1 * qJ(3,3);
t671 = -0.2e1 * t304;
t670 = -0.2e1 * t305;
t669 = -0.2e1 * t306;
t668 = 0.2e1 * t321;
t667 = 0.2e1 * t323;
t666 = 0.2e1 * t325;
t665 = -0.2e1 * t331;
t664 = pkin(2) * g(1);
t663 = pkin(2) * g(2);
t662 = pkin(5) * g(3);
t660 = t46 * pkin(2);
t659 = t47 * pkin(2);
t657 = t48 * pkin(2);
t655 = t671 + 0.1e1;
t654 = t670 + 0.1e1;
t653 = t669 + 0.1e1;
t649 = pkin(2) * t316;
t648 = pkin(2) * t318;
t647 = pkin(2) * t320;
t637 = t46 * qJ(3,3);
t636 = t47 * qJ(3,2);
t635 = t48 * qJ(3,1);
t291 = g(3) * t322;
t199 = t202 ^ 2;
t613 = t199 * t230;
t97 = pkin(5) * t100;
t88 = pkin(1) * t613 + t291 - t97;
t85 = t88 + t695;
t634 = t85 * t315;
t292 = g(3) * t324;
t200 = t203 ^ 2;
t611 = t200 * t233;
t98 = pkin(5) * t101;
t89 = pkin(1) * t611 + t292 - t98;
t86 = t89 + t694;
t633 = t86 * t317;
t293 = g(3) * t326;
t201 = t204 ^ 2;
t609 = t201 * t236;
t99 = pkin(5) * t102;
t90 = pkin(1) * t609 + t293 - t99;
t87 = t90 + t693;
t632 = t87 * t319;
t631 = qJ(3,1) * t102;
t629 = qJ(3,2) * t101;
t627 = qJ(3,3) * t100;
t625 = t100 * t334;
t624 = t101 * t337;
t623 = t102 * t340;
t619 = t181 * t334;
t618 = t182 * t334;
t617 = t183 * t337;
t616 = t184 * t337;
t615 = t185 * t340;
t614 = t186 * t340;
t612 = t199 * t315;
t610 = t200 * t317;
t608 = t201 * t319;
t594 = t205 * t334;
t592 = t206 * t337;
t590 = t207 * t340;
t238 = t282 * g(1) + t285 * g(2);
t577 = t238 * t315;
t576 = t238 * t321;
t239 = t283 * g(1) + t286 * g(2);
t575 = t239 * t317;
t574 = t239 * t323;
t240 = t284 * g(1) + t287 * g(2);
t573 = t240 * t319;
t572 = t240 * t325;
t568 = t282 * t331;
t564 = t283 * t331;
t560 = t284 * t331;
t556 = t285 * t331;
t552 = t286 * t331;
t548 = t287 * t331;
t544 = t315 * t321;
t540 = t317 * t323;
t536 = t319 * t325;
t531 = t321 * t334;
t529 = t323 * t337;
t527 = t325 * t340;
t517 = t330 * t331;
t516 = t331 * t115;
t515 = t331 * t116;
t514 = t331 * t117;
t462 = t160 * t586;
t127 = 0.2e1 * t139 * t462;
t193 = pkin(2) * t613;
t501 = t193 + t127;
t460 = t161 * t582;
t128 = 0.2e1 * t140 * t460;
t194 = pkin(2) * t611;
t500 = t194 + t128;
t458 = t162 * t578;
t129 = 0.2e1 * t141 * t458;
t195 = pkin(2) * t609;
t499 = t195 + t129;
t344 = pkin(2) ^ 2;
t495 = -t333 + t344;
t494 = -t336 + t344;
t493 = -t339 + t344;
t491 = -0.4e1 * pkin(1) * t331;
t489 = qJ(3,1) * t665;
t488 = qJ(3,2) * t665;
t487 = qJ(3,3) * t665;
t484 = t46 * t588;
t483 = t47 * t584;
t482 = t48 * t580;
t478 = t100 * t543;
t477 = t100 * t531;
t476 = t101 * t539;
t475 = t101 * t529;
t474 = t102 * t535;
t473 = t102 * t527;
t472 = (0.4e1 * t154 - 0.2e1 * t136) * t607;
t471 = t122 * t607;
t470 = (0.4e1 * t155 - 0.2e1 * t137) * t603;
t469 = t124 * t603;
t468 = (0.4e1 * t156 - 0.2e1 * t138) * t599;
t467 = t126 * t599;
t457 = t230 * t612;
t456 = t233 * t610;
t455 = t236 * t608;
t453 = t160 * t606;
t449 = t161 * t602;
t445 = t162 * t598;
t442 = t205 * t531;
t441 = t206 * t529;
t440 = t207 * t527;
t421 = t100 * t544;
t419 = t101 * t540;
t417 = t102 * t536;
t211 = t291 + t695;
t212 = t292 + t694;
t213 = t293 + t693;
t403 = t199 * t437;
t402 = t200 * t434;
t401 = t201 * t431;
t175 = t655 * t613;
t176 = t654 * t611;
t177 = t653 * t609;
t399 = t193 * t671;
t398 = t194 * t670;
t397 = t195 * t669;
t390 = t321 * t403;
t389 = t323 * t402;
t388 = t325 * t401;
t109 = -qJ(3,3) * t160 * t229 + t516;
t151 = -0.2e1 * t154;
t214 = t315 * t415 + t542;
t217 = t415 + 0.2e1 * t420;
t169 = (-t282 * t429 + t285 * t487) * t304 + (-t217 * t568 + t285 * t680) * t321 - t214 * t571 + t270 * t556;
t170 = (t282 * t487 + t285 * t429) * t304 + (t217 * t556 + t282 * t680) * t321 + t214 * t559 + t270 * t568;
t187 = t304 * t428 + ((t701 + pkin(1)) * t322 + t276) * t510 + qJ(3,3) * (t270 * t322 + t315 * t276);
t308 = t331 ^ 2;
t372 = -((-t513 * t112 + (-t202 * (0.3e1 * t333 + t400) * t331 + (t202 * t491 - t622) * t279) * t229) * t321 - (t250 * t315 * t607 + t504) * qJ(3,3) + (t330 * (t136 + t151 - 0.2e1 * t409) * qJ(3,3) + (-(t308 - 0.3e1 * t333) * t510 - 0.3e1 * (-t333 / 0.3e1 + t308) * t279 + (t333 - t308) * t692) * t607) * t304) * t452 - ((t109 * t331 + t430 * t454) * t321 + pkin(1) * t516 + (t109 * t315 + (t655 * t202 * t517 - pkin(1) * t160) * t229) * qJ(3,3)) * t462 - (-t145 * t510 + ((pkin(1) * t334 + t315) * t331 * t160 + (t304 - 0.1e1) * qJ(3,3) * t605) * t229) * t139 * t586 + (t169 * t313 + t170 * t314 + t187 * t312) * t588;
t110 = -qJ(3,2) * t161 * t232 + t515;
t152 = -0.2e1 * t155;
t215 = t317 * t414 + t538;
t219 = t414 + 0.2e1 * t418;
t171 = (-t283 * t426 + t286 * t488) * t305 + (-t219 * t564 + t286 * t681) * t323 - t215 * t567 + t271 * t552;
t172 = (t283 * t488 + t286 * t426) * t305 + (t219 * t552 + t283 * t681) * t323 + t215 * t555 + t271 * t564;
t188 = t305 * t425 + ((t700 + pkin(1)) * t324 + t277) * t509 + qJ(3,2) * (t271 * t324 + t317 * t277);
t371 = -((-t512 * t113 + (-t203 * (0.3e1 * t336 + t400) * t331 + (t203 * t491 - t621) * t280) * t232) * t323 - (t251 * t317 * t603 + t503) * qJ(3,2) + (t330 * (t137 + t152 - 0.2e1 * t408) * qJ(3,2) + (-(t308 - 0.3e1 * t336) * t509 - 0.3e1 * (-t336 / 0.3e1 + t308) * t280 + (t336 - t308) * t692) * t603) * t305) * t448 - ((t110 * t331 + t427 * t450) * t323 + pkin(1) * t515 + (t110 * t317 + (t654 * t203 * t517 - pkin(1) * t161) * t232) * qJ(3,2)) * t460 - (-t146 * t509 + ((pkin(1) * t337 + t317) * t331 * t161 + (t305 - 0.1e1) * qJ(3,2) * t601) * t232) * t140 * t582 + (t171 * t313 + t172 * t314 + t188 * t312) * t584;
t111 = -qJ(3,1) * t162 * t235 + t514;
t153 = -0.2e1 * t156;
t216 = t319 * t413 + t534;
t221 = t413 + 0.2e1 * t416;
t173 = (-t284 * t423 + t287 * t489) * t306 + (-t221 * t560 + t287 * t682) * t325 - t216 * t563 + t272 * t548;
t174 = (t284 * t489 + t287 * t423) * t306 + (t221 * t548 + t284 * t682) * t325 + t216 * t551 + t272 * t560;
t189 = t306 * t422 + ((t699 + pkin(1)) * t326 + t278) * t508 + qJ(3,1) * (t272 * t326 + t319 * t278);
t370 = -((-t511 * t114 + (-t204 * (0.3e1 * t339 + t400) * t331 + (t204 * t491 - t620) * t281) * t235) * t325 - (t252 * t319 * t599 + t502) * qJ(3,1) + (t330 * (t138 + t153 - 0.2e1 * t407) * qJ(3,1) + (-(t308 - 0.3e1 * t339) * t508 - 0.3e1 * (-t339 / 0.3e1 + t308) * t281 + (t339 - t308) * t692) * t599) * t306) * t444 - ((t111 * t331 + t424 * t446) * t325 + pkin(1) * t514 + (t111 * t319 + (t653 * t204 * t517 - pkin(1) * t162) * t235) * qJ(3,1)) * t458 - (-t147 * t508 + ((pkin(1) * t340 + t319) * t331 * t162 + (t306 - 0.1e1) * qJ(3,1) * t597) * t235) * t141 * t578 + (t173 * t313 + t174 * t314 + t189 * t312) * t580;
t366 = -t372 + t660;
t365 = -t371 + t659;
t364 = -t370 + t657;
t363 = t372 + t576;
t362 = t371 + t574;
t361 = t370 + t572;
t332 = pkin(1) * g(3);
t255 = -t319 * pkin(2) + t528;
t254 = -t317 * pkin(2) + t530;
t253 = -t315 * pkin(2) + t532;
t144 = (-t159 * t341 + t201 * t306 - t201) * t236;
t143 = (-t158 * t338 + t200 * t305 - t200) * t233;
t142 = (-t157 * t335 + t199 * t304 - t199) * t230;
t93 = 0.2e1 * t631;
t92 = 0.2e1 * t629;
t91 = 0.2e1 * t627;
t81 = t87 * t325 + t573;
t80 = -t572 + t632;
t79 = t86 * t323 + t575;
t78 = -t574 + t633;
t77 = t85 * t321 + t577;
t76 = -t576 + t634;
t69 = t102 * t678 + 0.4e1 * t445;
t68 = t101 * t678 + 0.4e1 * t449;
t67 = t100 * t678 + 0.4e1 * t453;
t66 = (t102 * t319 + t404 * t666) * t319;
t65 = (t101 * t317 + t405 * t667) * t317;
t64 = (t100 * t315 + t406 * t668) * t315;
t63 = 0.2e1 * t417 + (t689 - 0.2e1) * t404;
t62 = 0.2e1 * t419 + (t690 - 0.2e1) * t405;
t61 = 0.2e1 * t421 + (t691 - 0.2e1) * t406;
t45 = -t358 * t580 + (t385 * t341 + (-t105 * t596 - t608) * t325) * t236;
t44 = -t359 * t584 + (t386 * t338 + (-t104 * t600 - t610) * t323) * t233;
t43 = -t360 * t588 + (t387 * t335 + (-t103 * t604 - t612) * t321) * t230;
t42 = -t319 * t464 + t48 * t325;
t41 = t48 * t319 + t325 * t464;
t40 = -t317 * t465 + t47 * t323;
t39 = t47 * t317 + t323 * t465;
t38 = -t315 * t466 + t46 * t321;
t37 = t46 * t315 + t321 * t466;
t30 = t325 * t696 - t705;
t29 = -t36 * t325 - t704;
t28 = t323 * t697 - t706;
t27 = -t35 * t323 - t703;
t26 = t321 * t698 - t707;
t25 = -t321 * t34 - t702;
t24 = t397 + (t455 * t674 - t87) * t325 - t573 + 0.2e1 * t635 + t499;
t23 = t398 + (t456 * t673 - t86) * t323 - t575 + 0.2e1 * t636 + t500;
t22 = t399 + (t457 * t672 - t85) * t321 - t577 + 0.2e1 * t637 + t501;
t21 = t69 * t306 + ((t93 - t468) * t319 + t696) * t325 - t705 - 0.2e1 * t445;
t20 = t68 * t305 + ((t92 - t470) * t317 + t697) * t323 - t706 - 0.2e1 * t449;
t19 = t67 * t304 + ((t91 - t472) * t315 + t698) * t321 - t707 - 0.2e1 * t453;
t18 = (t468 - 0.2e1 * t631) * t306 + (t69 * t319 + t36) * t325 + t704 + (t153 + 0.2e1 * t138) * t599 + t93;
t17 = (t470 - 0.2e1 * t629) * t305 + (t68 * t317 + t35) * t323 + t703 + (t152 + 0.2e1 * t137) * t603 + t92;
t16 = (t472 - 0.2e1 * t627) * t304 + (t67 * t315 + t34) * t321 + t702 + (t151 + 0.2e1 * t136) * t607 + t91;
t15 = (t195 * t666 + t87) * t319 + 0.2e1 * t657 + qJ(3,1) * t177 - t361;
t14 = (t194 * t667 + t86) * t317 + 0.2e1 * t659 + qJ(3,2) * t176 - t362;
t13 = (t193 * t668 + t85) * t315 + 0.2e1 * t660 + qJ(3,3) * t175 - t363;
t12 = -t632 - t657 + (-t159 * t340 + (-pkin(2) * t536 - qJ(3,1) + t630) * t201) * t236 + t361;
t11 = -t633 - t659 + (-t158 * t337 + (-pkin(2) * t540 - qJ(3,2) + t628) * t200) * t233 + t362;
t10 = -t634 - t660 + (-t157 * t334 + (-pkin(2) * t544 - qJ(3,3) + t626) * t199) * t230 + t363;
t9 = (-pkin(2) * t464 + t129 + t635) * t325 + (-t159 * t579 - t364) * t319 + 0.2e1 * t99 - t213;
t8 = (-pkin(2) * t465 + t128 + t636) * t323 + (-t158 * t583 - t365) * t317 + 0.2e1 * t98 - t212;
t7 = (-pkin(2) * t466 + t127 + t637) * t321 + (-t157 * t587 - t366) * t315 + 0.2e1 * t97 - t211;
t6 = qJ(3,1) * t397 + ((-g(1) * t534 - t663) * t287 + (g(2) * t534 - t664) * t284 - qJ(3,1) * t90 + t493 * t455) * t325 + ((g(1) * t647 - g(2) * qJ(3,1)) * t287 + (-g(1) * qJ(3,1) - g(2) * t647) * t284 + pkin(2) * t90) * t319 + t339 * t48 + t499 * qJ(3,1) + pkin(2) * t364;
t5 = qJ(3,2) * t398 + ((-g(1) * t538 - t663) * t286 + (g(2) * t538 - t664) * t283 - qJ(3,2) * t89 + t494 * t456) * t323 + ((g(1) * t648 - g(2) * qJ(3,2)) * t286 + (-g(1) * qJ(3,2) - g(2) * t648) * t283 + pkin(2) * t89) * t317 + t336 * t47 + t500 * qJ(3,2) + pkin(2) * t365;
t4 = qJ(3,3) * t399 + ((-g(1) * t542 - t663) * t285 + (g(2) * t542 - t664) * t282 - qJ(3,3) * t88 + t495 * t457) * t321 + ((g(1) * t649 - g(2) * qJ(3,3)) * t285 + (-g(1) * qJ(3,3) - g(2) * t649) * t282 + pkin(2) * t88) * t315 + t333 * t46 + t501 * qJ(3,3) + pkin(2) * t366;
t3 = (t156 - t138 / 0.2e1) * t481 * t689 + (pkin(5) * t129 + (t445 * t686 + t36) * qJ(3,1) + ((-t467 + t631) * t686 + t696) * pkin(2)) * t325 + (t141 * t444 * t692 - t36 * pkin(2) + pkin(5) * t370 + qJ(3,1) * t696) * t319 + (-pkin(1) * t243 - t662) * t326 + (-pkin(5) * t243 + t332) * t320 + t467 * t674 + (t493 * t306 + t297 + t339) * t102;
t2 = (t155 - t137 / 0.2e1) * t480 * t690 + (pkin(5) * t128 + (t449 * t687 + t35) * qJ(3,2) + ((-t469 + t629) * t687 + t697) * pkin(2)) * t323 + (t140 * t448 * t692 - t35 * pkin(2) + pkin(5) * t371 + qJ(3,2) * t697) * t317 + (-pkin(1) * t242 - t662) * t324 + (-pkin(5) * t242 + t332) * t318 + t469 * t673 + (t494 * t305 + t297 + t336) * t101;
t1 = (t154 - t136 / 0.2e1) * t479 * t691 + (pkin(5) * t127 + (t453 * t688 + t34) * qJ(3,3) + ((-t471 + t627) * t688 + t698) * pkin(2)) * t321 + (t139 * t452 * t692 - t34 * pkin(2) + pkin(5) * t372 + qJ(3,3) * t698) * t315 + (-pkin(1) * t241 - t662) * t322 + (-pkin(5) * t241 + t332) * t316 + t471 * t672 + (t495 * t304 + t297 + t333) * t100;
t31 = [t100 * t438 + t101 * t435 + t102 * t432, t208 * t438 + t209 * t435 + t210 * t432, t211 * t438 + t212 * t435 + t213 * t432, -t182 * t390 - t184 * t389 - t186 * t388 + t432 * t66 + t435 * t65 + t438 * t64, (t177 * t614 + t550 * t63) * t235 + (t176 * t616 + t554 * t62) * t232 + (t175 * t618 + t558 * t61) * t229, (t186 * t474 + t41 * t550) * t235 + (t184 * t476 + t39 * t554) * t232 + (t182 * t478 + t37 * t558) * t229, (t186 * t473 + t42 * t550) * t235 + (t184 * t475 + t40 * t554) * t232 + (t182 * t477 + t38 * t558) * t229, t182 * t484 + t184 * t483 + t186 * t482, (t30 * t550 + t614 * t80) * t235 + (t28 * t554 + t616 * t78) * t232 + (t26 * t558 + t618 * t76) * t229, (t29 * t550 + t614 * t81) * t235 + (t27 * t554 + t616 * t79) * t232 + (t25 * t558 + t618 * t77) * t229, (t21 * t550 + (t15 * t186 + t174 * t45) * t340) * t235 + (t20 * t554 + (t14 * t184 + t172 * t44) * t337) * t232 + (t19 * t558 + (t13 * t182 + t170 * t43) * t334) * t229, (t9 * t550 + (t174 * t319 + t186 * t255) * t623) * t235 + (t8 * t554 + (t172 * t317 + t184 * t254) * t624) * t232 + (t7 * t558 + (t170 * t315 + t182 * t253) * t625) * t229, (t18 * t550 + (t144 * t174 + t186 * t24) * t340) * t235 + (t17 * t554 + (t143 * t172 + t184 * t23) * t337) * t232 + (t16 * t558 + (t142 * t170 + t182 * t22) * t334) * t229, (t3 * t550 + (t12 * t174 + t186 * t6) * t340) * t235 + (t2 * t554 + (t11 * t172 + t184 * t5) * t337) * t232 + (t1 * t558 + (t10 * t170 + t182 * t4) * t334) * t229, t314 - g(1); -t100 * t439 - t101 * t436 - t102 * t433, -t208 * t439 - t209 * t436 - t210 * t433, -t211 * t439 - t212 * t436 - t213 * t433, -t181 * t390 - t183 * t389 - t185 * t388 - t433 * t66 - t436 * t65 - t439 * t64, (t177 * t615 - t562 * t63) * t235 + (t176 * t617 - t566 * t62) * t232 + (t175 * t619 - t570 * t61) * t229, (t185 * t474 - t41 * t562) * t235 + (t183 * t476 - t39 * t566) * t232 + (t181 * t478 - t37 * t570) * t229, (t185 * t473 - t42 * t562) * t235 + (t183 * t475 - t40 * t566) * t232 + (t181 * t477 - t38 * t570) * t229, t181 * t484 + t183 * t483 + t185 * t482, (-t30 * t562 + t615 * t80) * t235 + (-t28 * t566 + t617 * t78) * t232 + (-t26 * t570 + t619 * t76) * t229, (-t29 * t562 + t615 * t81) * t235 + (-t27 * t566 + t617 * t79) * t232 + (-t25 * t570 + t619 * t77) * t229, (-t21 * t562 + (t15 * t185 + t173 * t45) * t340) * t235 + (-t20 * t566 + (t14 * t183 + t171 * t44) * t337) * t232 + (-t19 * t570 + (t13 * t181 + t169 * t43) * t334) * t229, (-t9 * t562 + (t173 * t319 + t185 * t255) * t623) * t235 + (-t8 * t566 + (t171 * t317 + t183 * t254) * t624) * t232 + (-t7 * t570 + (t169 * t315 + t181 * t253) * t625) * t229, (-t18 * t562 + (t144 * t173 + t185 * t24) * t340) * t235 + (-t17 * t566 + (t143 * t171 + t183 * t23) * t337) * t232 + (-t16 * t570 + (t142 * t169 + t181 * t22) * t334) * t229, (-t3 * t562 + (t12 * t173 + t185 * t6) * t340) * t235 + (-t2 * t566 + (t11 * t171 + t183 * t5) * t337) * t232 + (-t1 * t570 + (t10 * t169 + t181 * t4) * t334) * t229, t313 - g(2); -t100 * t589 - t101 * t585 - t102 * t581, -t208 * t589 - t209 * t585 - t210 * t581, -t211 * t589 - t212 * t585 - t213 * t581, -t205 * t304 * t403 - t206 * t305 * t402 - t207 * t306 * t401 - t581 * t66 - t585 * t65 - t589 * t64, (t177 * t440 - t320 * t63) * t235 + (t176 * t441 - t318 * t62) * t232 + (t175 * t442 - t316 * t61) * t229, (-t320 * t41 + t417 * t590) * t235 + (-t318 * t39 + t419 * t592) * t232 + (-t316 * t37 + t421 * t594) * t229, (t102 * t306 * t590 - t320 * t42) * t235 + (t101 * t305 * t592 - t318 * t40) * t232 + (t100 * t304 * t594 - t316 * t38) * t229, t229 * t442 * t46 + t232 * t441 * t47 + t235 * t440 * t48, (-t30 * t320 + t440 * t80) * t235 + (-t28 * t318 + t441 * t78) * t232 + (-t26 * t316 + t442 * t76) * t229, (-t29 * t320 + t440 * t81) * t235 + (-t27 * t318 + t441 * t79) * t232 + (-t25 * t316 + t442 * t77) * t229, (-t21 * t320 + (t15 * t591 + t189 * t45) * t340) * t235 + (-t20 * t318 + (t14 * t593 + t188 * t44) * t337) * t232 + (-t19 * t316 + (t13 * t595 + t187 * t43) * t334) * t229, (-t320 * t9 + (t189 * t319 + t255 * t591) * t623) * t235 + (-t318 * t8 + (t188 * t317 + t254 * t593) * t624) * t232 + (-t316 * t7 + (t187 * t315 + t253 * t595) * t625) * t229, (-t18 * t320 + (t144 * t189 + t24 * t591) * t340) * t235 + (-t17 * t318 + (t143 * t188 + t23 * t593) * t337) * t232 + (-t16 * t316 + (t142 * t187 + t22 * t595) * t334) * t229, (-t3 * t320 + (t12 * t189 + t591 * t6) * t340) * t235 + (-t2 * t318 + (t11 * t188 + t5 * t593) * t337) * t232 + (-t1 * t316 + (t10 * t187 + t4 * t595) * t334) * t229, t312 - g(3);];
tauX_reg  = t31;