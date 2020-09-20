% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR12V1G3A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRPRR12V1G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:10:12
% EndTime: 2020-08-06 19:10:29
% DurationCPUTime: 16.77s
% Computational Cost: add. (71811->663), mult. (123306->1196), div. (7551->12), fcn. (78729->18), ass. (0->492)
t241 = legFrame(3,2);
t225 = cos(t241);
t531 = t225 * g(1);
t222 = sin(t241);
t534 = t222 * g(2);
t181 = t531 - t534;
t248 = sin(qJ(1,3));
t231 = g(3) * t248;
t254 = cos(qJ(1,3));
t586 = t181 * t254 - t231;
t242 = legFrame(2,2);
t226 = cos(t242);
t530 = t226 * g(1);
t223 = sin(t242);
t533 = t223 * g(2);
t182 = t530 - t533;
t250 = sin(qJ(1,2));
t232 = g(3) * t250;
t256 = cos(qJ(1,2));
t585 = t182 * t256 - t232;
t243 = legFrame(1,2);
t227 = cos(t243);
t529 = t227 * g(1);
t224 = sin(t243);
t532 = t224 * g(2);
t183 = t529 - t532;
t252 = sin(qJ(1,1));
t233 = g(3) * t252;
t258 = cos(qJ(1,1));
t584 = t183 * t258 - t233;
t260 = xDP(2);
t262 = pkin(1) + pkin(2);
t220 = t260 * t262;
t261 = xDP(1);
t221 = t261 * t262;
t259 = xDP(3);
t399 = t259 * t262;
t508 = qJ(3,1) * t261;
t509 = qJ(3,1) * t260;
t120 = (t221 * t258 - t509) * t227 + (-t220 * t258 - t508) * t224 - t252 * t399;
t526 = t252 * pkin(4);
t219 = t260 * t526;
t257 = cos(qJ(2,1));
t239 = t257 ^ 2;
t251 = sin(qJ(2,1));
t384 = t261 * t526;
t230 = t258 * pkin(4);
t415 = t251 * t252;
t188 = qJ(3,1) * t415 + t230;
t400 = t259 * t188;
t563 = -t258 * t509 + t221;
t566 = t258 * t508 + t220;
t102 = t120 * t239 + ((t566 * t251 - t384) * t227 + (t563 * t251 + t219) * t224 - t400) * t257 + qJ(3,1) * (t224 * t261 + t227 * t260);
t413 = t251 * t258;
t189 = qJ(3,1) * t413 - t526;
t396 = t262 * t251;
t403 = t258 * t262;
t431 = t224 * qJ(3,1);
t114 = (t227 * t403 - t431) * t239 + (t189 * t227 + t224 * t396) * t257 + t431;
t425 = t227 * qJ(3,1);
t117 = (-t224 * t403 - t425) * t239 + (-t224 * t189 + t227 * t396) * t257 + t425;
t393 = t262 * t257;
t165 = t252 * t393 + t188;
t204 = t251 * qJ(3,1) + t393;
t196 = 0.1e1 / t204;
t244 = xDDP(3);
t245 = xDDP(2);
t246 = xDDP(1);
t270 = 0.1e1 / qJ(3,1);
t439 = t196 * t251;
t108 = -t227 * t384 + t120 * t257 - t400 + t219 * t224 + (t563 * t224 + t566 * t227) * t251;
t271 = 0.1e1 / qJ(3,1) ^ 2;
t493 = t108 * t271;
t150 = (t224 * t260 - t227 * t261) * t252 - t258 * t259;
t535 = pkin(4) * t150;
t135 = t439 * t535;
t437 = t196 * t270;
t346 = t102 * t437;
t90 = t262 * t346 + t135;
t580 = ((t165 * t244 * t257 - t114 * t246 - t117 * t245) * t270 + (t102 * t439 + t90 * t257) * t493) * t196;
t505 = qJ(3,2) * t261;
t506 = qJ(3,2) * t260;
t119 = (t221 * t256 - t506) * t226 + (-t220 * t256 - t505) * t223 - t250 * t399;
t527 = t250 * pkin(4);
t218 = t260 * t527;
t255 = cos(qJ(2,2));
t238 = t255 ^ 2;
t249 = sin(qJ(2,2));
t385 = t261 * t527;
t229 = t256 * pkin(4);
t419 = t249 * t250;
t186 = qJ(3,2) * t419 + t229;
t401 = t259 * t186;
t564 = -t256 * t506 + t221;
t567 = t256 * t505 + t220;
t101 = t119 * t238 + ((t567 * t249 - t385) * t226 + (t564 * t249 + t218) * t223 - t401) * t255 + qJ(3,2) * (t223 * t261 + t226 * t260);
t417 = t249 * t256;
t187 = qJ(3,2) * t417 - t527;
t397 = t262 * t249;
t406 = t256 * t262;
t433 = t223 * qJ(3,2);
t113 = (t226 * t406 - t433) * t238 + (t187 * t226 + t223 * t397) * t255 + t433;
t427 = t226 * qJ(3,2);
t116 = (-t223 * t406 - t427) * t238 + (-t223 * t187 + t226 * t397) * t255 + t427;
t394 = t262 * t255;
t164 = t250 * t394 + t186;
t203 = t249 * qJ(3,2) + t394;
t193 = 0.1e1 / t203;
t267 = 0.1e1 / qJ(3,2);
t443 = t193 * t249;
t107 = -t226 * t385 + t119 * t255 - t401 + t218 * t223 + (t564 * t223 + t567 * t226) * t249;
t268 = 0.1e1 / qJ(3,2) ^ 2;
t494 = t107 * t268;
t149 = (t223 * t260 - t226 * t261) * t250 - t256 * t259;
t536 = pkin(4) * t149;
t134 = t443 * t536;
t441 = t193 * t267;
t348 = t101 * t441;
t89 = t262 * t348 + t134;
t579 = (t267 * (t164 * t244 * t255 - t113 * t246 - t116 * t245) + (t101 * t443 + t89 * t255) * t494) * t193;
t502 = qJ(3,3) * t261;
t503 = qJ(3,3) * t260;
t118 = (t221 * t254 - t503) * t225 + (-t220 * t254 - t502) * t222 - t248 * t399;
t528 = t248 * pkin(4);
t217 = t260 * t528;
t253 = cos(qJ(2,3));
t237 = t253 ^ 2;
t247 = sin(qJ(2,3));
t386 = t261 * t528;
t228 = t254 * pkin(4);
t423 = t247 * t248;
t184 = qJ(3,3) * t423 + t228;
t402 = t259 * t184;
t565 = -t254 * t503 + t221;
t568 = t254 * t502 + t220;
t100 = t118 * t237 + ((t568 * t247 - t386) * t225 + (t565 * t247 + t217) * t222 - t402) * t253 + qJ(3,3) * (t222 * t261 + t225 * t260);
t421 = t247 * t254;
t185 = qJ(3,3) * t421 - t528;
t398 = t262 * t247;
t409 = t254 * t262;
t435 = t222 * qJ(3,3);
t112 = (t225 * t409 - t435) * t237 + (t185 * t225 + t222 * t398) * t253 + t435;
t429 = t225 * qJ(3,3);
t115 = (-t222 * t409 - t429) * t237 + (-t222 * t185 + t225 * t398) * t253 + t429;
t395 = t262 * t253;
t163 = t248 * t395 + t184;
t202 = t247 * qJ(3,3) + t395;
t190 = 0.1e1 / t202;
t264 = 0.1e1 / qJ(3,3);
t447 = t190 * t247;
t106 = -t225 * t386 + t118 * t253 - t402 + t217 * t222 + (t565 * t222 + t568 * t225) * t247;
t265 = 0.1e1 / qJ(3,3) ^ 2;
t495 = t106 * t265;
t148 = (t222 * t260 - t225 * t261) * t248 - t254 * t259;
t537 = pkin(4) * t148;
t133 = t447 * t537;
t445 = t190 * t264;
t350 = t100 * t445;
t88 = t262 * t350 + t133;
t578 = ((t163 * t244 * t253 - t112 * t246 - t115 * t245) * t264 + (t100 * t447 + t88 * t253) * t495) * t190;
t577 = t202 * t254 - t528;
t576 = t203 * t256 - t527;
t575 = t204 * t258 - t526;
t194 = 0.1e1 / t203 ^ 2;
t440 = t194 * t267;
t330 = t149 * t440;
t306 = t101 * t330;
t418 = t249 * t255;
t104 = t107 * t267;
t195 = t193 * t194;
t408 = t255 * qJ(3,2);
t291 = -t397 + t408;
t426 = t226 * t250;
t314 = t193 * t426;
t432 = t223 * t250;
t315 = t193 * t432;
t391 = t262 * t267;
t442 = t193 * t256;
t463 = t149 * t267;
t464 = t149 * t194;
t50 = -t246 * t314 + t245 * t315 - t244 * t442 - (t104 * t249 + (-t536 + (-t249 * t391 + t255) * t101) * t193) * t464 - t291 * t195 * t101 * t463 - t249 * t107 * t330;
t542 = 0.2e1 * t238 - 0.1e1;
t37 = t542 * t306 + t50 * t418;
t574 = -0.2e1 * t37;
t191 = 0.1e1 / t202 ^ 2;
t444 = t191 * t264;
t333 = t148 * t444;
t307 = t100 * t333;
t422 = t247 * t253;
t103 = t106 * t264;
t192 = t190 * t191;
t411 = t253 * qJ(3,3);
t290 = -t398 + t411;
t428 = t225 * t248;
t316 = t190 * t428;
t434 = t222 * t248;
t317 = t190 * t434;
t392 = t262 * t264;
t446 = t190 * t254;
t466 = t148 * t264;
t467 = t148 * t191;
t49 = -t246 * t316 + t245 * t317 - t244 * t446 - (t103 * t247 + (-t537 + (-t247 * t392 + t253) * t100) * t190) * t467 - t290 * t192 * t100 * t466 - t247 * t106 * t333;
t543 = 0.2e1 * t237 - 0.1e1;
t38 = t543 * t307 + t49 * t422;
t573 = -0.2e1 * t38;
t197 = 0.1e1 / t204 ^ 2;
t436 = t197 * t270;
t327 = t150 * t436;
t305 = t102 * t327;
t414 = t251 * t257;
t105 = t108 * t270;
t198 = t196 * t197;
t405 = t257 * qJ(3,1);
t292 = -t396 + t405;
t424 = t227 * t252;
t312 = t196 * t424;
t430 = t224 * t252;
t313 = t196 * t430;
t390 = t262 * t270;
t438 = t196 * t258;
t460 = t150 * t270;
t461 = t150 * t197;
t51 = -t246 * t312 + t245 * t313 - t244 * t438 - (t105 * t251 + (-t535 + (-t251 * t390 + t257) * t102) * t196) * t461 - t292 * t198 * t102 * t460 - t251 * t108 * t327;
t541 = 0.2e1 * t239 - 0.1e1;
t39 = t541 * t305 + t51 * t414;
t572 = -0.2e1 * t39;
t151 = t585 * t249;
t179 = t223 * g(1) + t226 * g(2);
t127 = -t179 * t255 + t151;
t152 = t586 * t247;
t178 = t222 * g(1) + t225 * g(2);
t129 = -t178 * t253 + t152;
t153 = t584 * t251;
t180 = t224 * g(1) + t227 * g(2);
t131 = -t180 * t257 + t153;
t234 = g(3) * t254;
t571 = t181 * t248 + t234;
t235 = g(3) * t256;
t570 = t182 * t250 + t235;
t236 = g(3) * t258;
t569 = t183 * t252 + t236;
t559 = 2 * pkin(1);
t558 = -0.2e1 * qJ(3,1);
t557 = -0.2e1 * qJ(3,2);
t556 = -0.2e1 * qJ(3,3);
t555 = 0.2e1 * t253;
t554 = 0.2e1 * t255;
t553 = 0.2e1 * t257;
t552 = -pkin(4) / 0.2e1;
t551 = pkin(1) * g(1);
t550 = pkin(1) * g(2);
t549 = pkin(1) * t49;
t548 = pkin(1) * t50;
t547 = pkin(1) * t51;
t263 = qJ(3,3) ^ 2;
t272 = pkin(4) ^ 2;
t468 = t148 * t190;
t335 = t237 * t468;
t501 = t100 * t264;
t92 = pkin(1) * t350;
t74 = t92 - t103;
t61 = pkin(2) * t350 + t74;
t522 = t247 * t61;
t375 = ((qJ(3,3) + t262) * (-qJ(3,3) + t262) * t335 + 0.2e1 * (t148 * t398 + t501 * t552) * t190 * t411 + pkin(4) * t522 + (t263 + t272) * t468) * t466;
t500 = t100 * t265;
t383 = (-(t133 + t61) * t395 + (pkin(4) * t335 - t522) * qJ(3,3)) * t500;
t34 = (t253 * t375 - t383) * t191 - t578;
t546 = t34 * pkin(1);
t266 = qJ(3,2) ^ 2;
t465 = t149 * t193;
t332 = t238 * t465;
t499 = t101 * t267;
t94 = pkin(1) * t348;
t76 = t94 - t104;
t62 = pkin(2) * t348 + t76;
t521 = t249 * t62;
t372 = ((qJ(3,2) + t262) * (-qJ(3,2) + t262) * t332 + 0.2e1 * (t149 * t397 + t499 * t552) * t193 * t408 + pkin(4) * t521 + (t266 + t272) * t465) * t463;
t498 = t101 * t268;
t382 = (-(t134 + t62) * t394 + (pkin(4) * t332 - t521) * qJ(3,2)) * t498;
t35 = (t255 * t372 - t382) * t194 - t579;
t545 = t35 * pkin(1);
t269 = qJ(3,1) ^ 2;
t462 = t150 * t196;
t329 = t239 * t462;
t497 = t102 * t270;
t96 = pkin(1) * t346;
t78 = t96 - t105;
t63 = pkin(2) * t346 + t78;
t520 = t251 * t63;
t369 = ((qJ(3,1) + t262) * (-qJ(3,1) + t262) * t329 + 0.2e1 * (t150 * t396 + t497 * t552) * t196 * t405 + pkin(4) * t520 + (t269 + t272) * t462) * t460;
t496 = t102 * t271;
t381 = (-(t135 + t63) * t393 + (pkin(4) * t329 - t520) * qJ(3,1)) * t496;
t36 = (t257 * t369 - t381) * t197 - t580;
t544 = t36 * pkin(1);
t540 = pkin(1) * t254;
t539 = pkin(1) * t256;
t538 = pkin(1) * t258;
t525 = qJ(3,1) * t51;
t524 = qJ(3,2) * t50;
t523 = qJ(3,3) * t49;
t97 = t100 ^ 2;
t519 = t265 * t97;
t98 = t101 ^ 2;
t518 = t268 * t98;
t99 = t102 ^ 2;
t517 = t271 * t99;
t516 = t34 * qJ(3,3);
t515 = t35 * qJ(3,2);
t514 = t36 * qJ(3,1);
t145 = t148 ^ 2;
t474 = t145 * t191;
t142 = pkin(1) * t474;
t349 = t190 * t500;
t79 = 0.2e1 * t106 * t349;
t513 = t142 + t79;
t146 = t149 ^ 2;
t472 = t146 * t194;
t143 = pkin(1) * t472;
t347 = t193 * t498;
t80 = 0.2e1 * t107 * t347;
t512 = t143 + t80;
t147 = t150 ^ 2;
t470 = t147 * t197;
t144 = pkin(1) * t470;
t345 = t196 * t496;
t81 = 0.2e1 * t108 * t345;
t511 = t144 + t81;
t510 = qJ(3,1) * t258;
t507 = qJ(3,2) * t256;
t504 = qJ(3,3) * t254;
t492 = t112 * t264;
t491 = t113 * t267;
t490 = t114 * t270;
t489 = t115 * t264;
t488 = t116 * t267;
t487 = t117 * t270;
t486 = t127 * t267;
t451 = t179 * t249;
t128 = t255 * t585 + t451;
t485 = t128 * t267;
t484 = t129 * t264;
t453 = t178 * t247;
t130 = t253 * t586 + t453;
t483 = t130 * t264;
t482 = t131 * t270;
t449 = t180 * t251;
t132 = t257 * t584 + t449;
t481 = t132 * t270;
t136 = -t290 * t222 + t225 * t577;
t480 = t136 * t264;
t137 = -t291 * t223 + t226 * t576;
t479 = t137 * t267;
t138 = -t292 * t224 + t227 * t575;
t478 = t138 * t270;
t139 = -t222 * t577 - t290 * t225;
t477 = t139 * t264;
t140 = -t223 * t576 - t291 * t226;
t476 = t140 * t267;
t141 = -t224 * t575 - t292 * t227;
t475 = t141 * t270;
t473 = t145 * t247;
t471 = t146 * t249;
t469 = t147 * t251;
t459 = t163 * t264;
t458 = t164 * t267;
t457 = t165 * t270;
t172 = -t202 * t248 - t228;
t456 = t172 * t264;
t173 = -t203 * t250 - t229;
t455 = t173 * t267;
t174 = -t204 * t252 - t230;
t454 = t174 * t270;
t420 = t247 * t264;
t416 = t249 * t267;
t412 = t251 * t270;
t410 = t253 * t264;
t407 = t255 * t267;
t404 = t257 * t270;
t274 = pkin(1) ^ 2;
t389 = -t263 + t274;
t388 = -t266 + t274;
t387 = -t269 + t274;
t377 = (0.4e1 * t92 - 0.2e1 * t103) * t468;
t376 = t74 * t468;
t374 = (0.4e1 * t94 - 0.2e1 * t104) * t465;
t373 = t76 * t465;
t371 = (0.4e1 * t96 - 0.2e1 * t105) * t462;
t370 = t78 * t462;
t208 = -t247 * pkin(1) + t411;
t368 = t208 * t264 * t49;
t209 = -t249 * pkin(1) + t408;
t367 = t209 * t267 * t50;
t210 = -t251 * pkin(1) + t405;
t366 = t210 * t270 * t51;
t365 = t49 * t420;
t364 = t50 * t416;
t363 = t51 * t412;
t362 = t49 * t410;
t361 = t50 * t407;
t360 = t51 * t404;
t359 = t191 * t519;
t358 = t194 * t518;
t357 = t197 * t517;
t344 = t112 * t445;
t343 = t113 * t441;
t342 = t114 * t437;
t341 = t115 * t445;
t340 = t116 * t441;
t339 = t117 * t437;
t338 = t191 * t473;
t337 = t194 * t471;
t336 = t197 * t469;
t334 = t100 * t467;
t331 = t101 * t464;
t328 = t102 * t461;
t326 = t571 * t423;
t325 = t571 * t248 * t253;
t324 = t569 * t415;
t323 = t569 * t252 * t257;
t322 = t570 * t419;
t321 = t570 * t250 * t255;
t320 = t163 * t410;
t319 = t164 * t407;
t318 = t165 * t404;
t311 = -t274 + (-(2 * pkin(1)) - pkin(2)) * pkin(2);
t310 = t49 * t320;
t309 = t50 * t319;
t308 = t51 * t318;
t304 = t145 * t192 * t420;
t303 = t146 * t195 * t416;
t302 = t147 * t198 * t412;
t301 = t190 * t320;
t300 = t193 * t319;
t299 = t196 * t318;
t109 = t543 * t474;
t110 = t542 * t472;
t111 = t541 * t470;
t298 = -0.2e1 * t237 * t142;
t297 = -0.2e1 * t238 * t143;
t296 = -0.2e1 * t239 * t144;
t295 = t253 * t304;
t294 = t255 * t303;
t293 = t257 * t302;
t289 = -t190 * t375 - t246 * t480 - t245 * t477 - t244 * t456 + (t106 * t392 + ((-t263 + t311) * t501 + t290 * t537) * t190) * t349 + t88 * t495;
t288 = -t193 * t372 - t246 * t479 - t245 * t476 - t244 * t455 + (t107 * t391 + ((-t266 + t311) * t499 + t291 * t536) * t193) * t347 + t89 * t494;
t287 = -t196 * t369 - t246 * t478 - t245 * t475 - t244 * t454 + (t108 * t390 + ((-t269 + t311) * t497 + t292 * t535) * t196) * t345 + t90 * t493;
t286 = t289 + t546;
t285 = t288 + t545;
t284 = t287 + t544;
t84 = (t147 * t239 - t147 - t517) * t197;
t83 = (t146 * t238 - t146 - t518) * t194;
t82 = (t145 * t237 - t145 - t519) * t191;
t48 = 0.2e1 * t525;
t47 = 0.2e1 * t524;
t46 = 0.2e1 * t523;
t42 = (t51 * t251 + t305 * t553) * t251;
t41 = (t50 * t249 + t306 * t554) * t249;
t40 = (t49 * t247 + t307 * t555) * t247;
t33 = (0.4e1 * t328 + 0.2e1 * t547) * t239 + ((t48 - t371) * t251 + t569) * t257 - 0.2e1 * t328;
t32 = (0.4e1 * t331 + 0.2e1 * t548) * t238 + ((t47 - t374) * t249 + t570) * t255 - 0.2e1 * t331;
t31 = (0.4e1 * t334 + 0.2e1 * t549) * t237 + ((t46 - t377) * t247 + t571) * t253 - 0.2e1 * t334;
t30 = (t381 + (-t369 - t469) * t257) * t197 + t580;
t29 = (t382 + (-t372 - t471) * t255) * t194 + t579;
t28 = (t383 + (-t375 - t473) * t253) * t191 + t578;
t27 = -t251 * t357 + t36 * t257;
t26 = t36 * t251 + t257 * t357;
t25 = -t249 * t358 + t35 * t255;
t24 = t35 * t249 + t255 * t358;
t23 = -t247 * t359 + t34 * t253;
t22 = t34 * t247 + t253 * t359;
t21 = t296 + (t336 * t558 - t584) * t257 - t449 + 0.2e1 * t514 + t511;
t20 = t297 + (t337 * t557 - t585) * t255 - t451 + 0.2e1 * t515 + t512;
t19 = t298 + (t338 * t556 - t586) * t253 - t453 + 0.2e1 * t516 + t513;
t18 = (t371 - 0.2e1 * t525) * t239 - 0.2e1 * t370 + t48 + ((0.2e1 * t328 + t547) * t553 + t569) * t251;
t17 = (t374 - 0.2e1 * t524) * t238 - 0.2e1 * t373 + t47 + ((0.2e1 * t331 + t548) * t554 + t570) * t249;
t16 = (t377 - 0.2e1 * t523) * t237 - 0.2e1 * t376 + t46 + ((0.2e1 * t334 + t549) * t555 + t571) * t247;
t15 = (t336 * t559 - t180) * t257 + t153 + 0.2e1 * t544 - qJ(3,1) * t111 + t287;
t14 = (t337 * t559 - t179) * t255 + t151 + 0.2e1 * t545 - qJ(3,2) * t110 + t288;
t13 = (t338 * t559 - t178) * t253 + t152 + 0.2e1 * t546 - qJ(3,3) * t109 + t289;
t12 = (-t270 * t99 + (-pkin(1) * t414 + qJ(3,1) * t239 - qJ(3,1)) * t147) * t197 - t284 - t131;
t11 = (-t267 * t98 + (-pkin(1) * t418 + qJ(3,2) * t238 - qJ(3,2)) * t146) * t194 - t285 - t127;
t10 = (-t264 * t97 + (-pkin(1) * t422 + qJ(3,3) * t237 - qJ(3,3)) * t145) * t191 - t286 - t129;
t9 = (0.4e1 * (t96 - t105 / 0.2e1) * qJ(3,1) * t462 + t387 * t51) * t239 + (0.2e1 * (qJ(3,1) * t328 + (-t370 + t525) * pkin(1)) * t251 + t569 * pkin(1)) * t257 + (((-t529 / 0.2e1 + t532 / 0.2e1) * t252 - t236 / 0.2e1) * t251 - t525 / 0.2e1 + t370) * t558;
t8 = (0.4e1 * (t94 - t104 / 0.2e1) * qJ(3,2) * t465 + t388 * t50) * t238 + (0.2e1 * (qJ(3,2) * t331 + (-t373 + t524) * pkin(1)) * t249 + t570 * pkin(1)) * t255 + (((-t530 / 0.2e1 + t533 / 0.2e1) * t250 - t235 / 0.2e1) * t249 - t524 / 0.2e1 + t373) * t557;
t7 = (0.4e1 * (t92 - t103 / 0.2e1) * qJ(3,3) * t468 + t389 * t49) * t237 + (0.2e1 * (qJ(3,3) * t334 + (-t376 + t523) * pkin(1)) * t247 + t571 * pkin(1)) * t253 + (((-t531 / 0.2e1 + t534 / 0.2e1) * t248 - t234 / 0.2e1) * t247 - t523 / 0.2e1 + t376) * t556;
t6 = qJ(3,1) * t296 + ((-g(1) * t510 - t550) * t227 + (g(2) * t510 - t551) * t224 + qJ(3,1) * t233 + t387 * t336) * t257 + ((g(1) * t538 - g(2) * qJ(3,1)) * t227 + (-g(1) * qJ(3,1) - g(2) * t538) * t224 - pkin(1) * t233) * t251 + t36 * t269 + t511 * qJ(3,1) + t284 * pkin(1);
t5 = qJ(3,2) * t297 + ((-g(1) * t507 - t550) * t226 + (g(2) * t507 - t551) * t223 + qJ(3,2) * t232 + t388 * t337) * t255 + ((g(1) * t539 - g(2) * qJ(3,2)) * t226 + (-g(1) * qJ(3,2) - g(2) * t539) * t223 - pkin(1) * t232) * t249 + t35 * t266 + t512 * qJ(3,2) + t285 * pkin(1);
t4 = qJ(3,3) * t298 + ((-g(1) * t504 - t550) * t225 + (g(2) * t504 - t551) * t222 + qJ(3,3) * t231 + t389 * t338) * t253 + ((g(1) * t540 - g(2) * qJ(3,3)) * t225 + (-g(1) * qJ(3,3) - g(2) * t540) * t222 - pkin(1) * t231) * t247 + t34 * t263 + t513 * qJ(3,3) + t286 * pkin(1);
t3 = (-pkin(1) * t357 + t514 + t81) * t257 + (-t436 * t99 - t284) * t251 - t584;
t2 = (-pkin(1) * t358 + t515 + t80) * t255 + (-t440 * t98 - t285) * t249 - t585;
t1 = (-pkin(1) * t359 + t516 + t79) * t253 + (-t444 * t97 - t286) * t247 - t586;
t43 = [-t312 * t51 - t314 * t50 - t316 * t49, -t312 * t569 - t314 * t570 - t316 * t571, -t312 * t584 - t314 * t585 - t316 * t586, -t112 * t295 - t113 * t294 - t114 * t293 - t312 * t42 - t314 * t41 - t316 * t40, -t109 * t344 - t110 * t343 - t111 * t342 + t312 * t572 + t314 * t574 + t316 * t573, (t114 * t363 - t26 * t424) * t196 + (t113 * t364 - t24 * t426) * t193 + (t112 * t365 - t22 * t428) * t190, (t114 * t360 - t27 * t424) * t196 + (t113 * t361 - t25 * t426) * t193 + (t112 * t362 - t23 * t428) * t190, t34 * t344 + t342 * t36 + t343 * t35, (t114 * t482 - t227 * t323) * t196 + (t113 * t486 - t226 * t321) * t193 + (t112 * t484 - t225 * t325) * t190, (t114 * t481 + t227 * t324) * t196 + (t113 * t485 + t226 * t322) * t193 + (t112 * t483 + t225 * t326) * t190, t28 * t480 + t29 * t479 + t30 * t478 + (t15 * t490 - t33 * t424) * t196 + (t14 * t491 - t32 * t426) * t193 + (t13 * t492 - t31 * t428) * t190, t136 * t365 + t137 * t364 + t138 * t363 + (t114 * t366 - t3 * t424) * t196 + (t113 * t367 - t2 * t426) * t193 + (-t1 * t428 + t112 * t368) * t190, t82 * t480 + t83 * t479 + t84 * t478 + (-t18 * t424 + t21 * t490) * t196 + (-t17 * t426 + t20 * t491) * t193 + (-t16 * t428 + t19 * t492) * t190, t10 * t480 + t11 * t479 + t12 * t478 + (-t424 * t9 + t490 * t6) * t196 + (-t426 * t8 + t491 * t5) * t193 + (t4 * t492 - t428 * t7) * t190, t246 - g(1); t313 * t51 + t315 * t50 + t317 * t49, t313 * t569 + t315 * t570 + t317 * t571, t313 * t584 + t315 * t585 + t317 * t586, -t115 * t295 - t116 * t294 - t117 * t293 + t313 * t42 + t315 * t41 + t317 * t40, -t109 * t341 - t110 * t340 - t111 * t339 + 0.2e1 * t313 * t39 + 0.2e1 * t315 * t37 + 0.2e1 * t317 * t38, (t117 * t363 + t26 * t430) * t196 + (t116 * t364 + t24 * t432) * t193 + (t115 * t365 + t22 * t434) * t190, (t117 * t360 + t27 * t430) * t196 + (t116 * t361 + t25 * t432) * t193 + (t115 * t362 + t23 * t434) * t190, t339 * t36 + t34 * t341 + t340 * t35, (t117 * t482 + t224 * t323) * t196 + (t116 * t486 + t223 * t321) * t193 + (t115 * t484 + t222 * t325) * t190, (t117 * t481 - t224 * t324) * t196 + (t116 * t485 - t223 * t322) * t193 + (t115 * t483 - t222 * t326) * t190, t28 * t477 + t29 * t476 + t30 * t475 + (t15 * t487 + t33 * t430) * t196 + (t14 * t488 + t32 * t432) * t193 + (t13 * t489 + t31 * t434) * t190, t139 * t365 + t140 * t364 + t141 * t363 + (t117 * t366 + t3 * t430) * t196 + (t116 * t367 + t2 * t432) * t193 + (t1 * t434 + t115 * t368) * t190, t82 * t477 + t83 * t476 + t84 * t475 + (t18 * t430 + t21 * t487) * t196 + (t17 * t432 + t20 * t488) * t193 + (t16 * t434 + t19 * t489) * t190, t10 * t477 + t11 * t476 + t12 * t475 + (t430 * t9 + t487 * t6) * t196 + (t432 * t8 + t488 * t5) * t193 + (t4 * t489 + t434 * t7) * t190, t245 - g(2); -t438 * t51 - t442 * t50 - t446 * t49, -t438 * t569 - t442 * t570 - t446 * t571, -t438 * t584 - t442 * t585 - t446 * t586, t163 * t237 * t304 + t164 * t238 * t303 + t165 * t239 * t302 - t40 * t446 - t41 * t442 - t42 * t438, t109 * t301 + t110 * t300 + t111 * t299 + t438 * t572 + t442 * t574 + t446 * t573, (-t251 * t308 - t258 * t26) * t196 + (-t24 * t256 - t249 * t309) * t193 + (-t22 * t254 - t247 * t310) * t190, (-t239 * t457 * t51 - t258 * t27) * t196 + (-t238 * t458 * t50 - t25 * t256) * t193 + (-t237 * t459 * t49 - t23 * t254) * t190, -t299 * t36 - t300 * t35 - t301 * t34, (-t131 * t457 - t258 * t569) * t257 * t196 + (-t127 * t458 - t256 * t570) * t255 * t193 + (-t129 * t459 - t254 * t571) * t253 * t190, (-t132 * t318 + t413 * t569) * t196 + (-t128 * t319 + t417 * t570) * t193 + (-t130 * t320 + t421 * t571) * t190, t28 * t456 + t29 * t455 + t30 * t454 + (-t15 * t318 - t258 * t33) * t196 + (-t14 * t319 - t256 * t32) * t193 + (-t13 * t320 - t254 * t31) * t190, t172 * t365 + t173 * t364 + t174 * t363 + (-t210 * t308 - t258 * t3) * t196 + (-t2 * t256 - t209 * t309) * t193 + (-t1 * t254 - t208 * t310) * t190, t82 * t456 + t83 * t455 + t84 * t454 + (-t18 * t258 - t21 * t318) * t196 + (-t17 * t256 - t20 * t319) * t193 + (-t16 * t254 - t19 * t320) * t190, t10 * t456 + t11 * t455 + t12 * t454 + (-t258 * t9 - t318 * t6) * t196 + (-t256 * t8 - t319 * t5) * t193 + (-t254 * t7 - t320 * t4) * t190, t244 - g(3);];
tauX_reg  = t43;
