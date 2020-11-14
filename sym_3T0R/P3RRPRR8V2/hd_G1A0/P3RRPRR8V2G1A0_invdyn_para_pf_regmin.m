% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x13]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRPRR8V2G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:03:06
% EndTime: 2020-08-06 21:03:21
% DurationCPUTime: 14.43s
% Computational Cost: add. (88397->574), mult. (151242->940), div. (14517->18), fcn. (84351->70), ass. (0->415)
t265 = sin(pkin(7));
t277 = sin(qJ(2,3));
t420 = t265 * t277;
t380 = pkin(3) * t420;
t283 = cos(qJ(2,3));
t266 = cos(pkin(7));
t220 = t266 * pkin(3);
t548 = t220 + pkin(2);
t424 = t548 * t283;
t109 = 0.1e1 / (-t380 + t424);
t271 = -qJ(3,3) - pkin(5);
t242 = -pkin(6) + t271;
t314 = t242 ^ 2;
t486 = pkin(3) * t265;
t112 = t277 * t548 + t283 * t486;
t293 = 0.2e1 * pkin(7);
t302 = pkin(3) ^ 2;
t499 = t302 / 0.2e1;
t131 = cos(t293) * t499 + pkin(2) * (t220 + pkin(2) / 0.2e1);
t290 = xDP(2);
t240 = pkin(1) * t290;
t291 = xDP(1);
t157 = -t242 * t291 + t240;
t238 = t291 * pkin(1);
t160 = t242 * t290 + t238;
t268 = legFrame(3,3);
t217 = t268 + qJ(1,3);
t295 = 0.2e1 * qJ(2,3);
t255 = sin(t295);
t258 = cos(t295);
t383 = 0.2e1 * t424;
t516 = 0.2e1 * t131;
t384 = t291 * t516;
t385 = t290 * t516;
t393 = 0.2e1 * pkin(3);
t198 = pkin(2) * t220;
t304 = pkin(2) ^ 2;
t241 = t302 + t304;
t556 = 0.2e1 * t198 + t241;
t396 = t556 * t291;
t397 = t556 * t290;
t166 = pkin(2) * t265 + sin(t293) * pkin(3) / 0.2e1;
t426 = t166 * t291;
t427 = t166 * t290;
t490 = pkin(3) * t166;
t289 = xDP(3);
t512 = 2 * t289;
t55 = (t258 * t384 + t160 * t383 + (-t160 * t420 - t255 * t426) * t393 + t396) * cos(t217) + (t258 * t385 + t157 * t383 + (-t157 * t420 - t255 * t427) * t393 + t397) * sin(t217) + (pkin(1) * t112 + t131 * t255 + t258 * t490) * t512;
t221 = sin(t268);
t224 = cos(t268);
t278 = sin(qJ(1,3));
t284 = cos(qJ(1,3));
t124 = -t221 * t278 + t224 * t284;
t127 = t284 * t221 + t278 * t224;
t227 = 0.1e1 / t242;
t445 = t109 * t112;
t76 = (t124 * t291 + t127 * t290 + t289 * t445) * t227;
t564 = t109 / t314 * t55 * t76;
t279 = sin(qJ(2,2));
t419 = t265 * t279;
t379 = pkin(3) * t419;
t285 = cos(qJ(2,2));
t423 = t548 * t285;
t110 = 0.1e1 / (-t379 + t423);
t272 = -qJ(3,2) - pkin(5);
t243 = -pkin(6) + t272;
t315 = t243 ^ 2;
t113 = t279 * t548 + t285 * t486;
t158 = -t243 * t291 + t240;
t161 = t243 * t290 + t238;
t269 = legFrame(2,3);
t218 = t269 + qJ(1,2);
t297 = 0.2e1 * qJ(2,2);
t256 = sin(t297);
t259 = cos(t297);
t382 = 0.2e1 * t423;
t56 = (t259 * t384 + t161 * t382 + (-t161 * t419 - t256 * t426) * t393 + t396) * cos(t218) + (t259 * t385 + t158 * t382 + (-t158 * t419 - t256 * t427) * t393 + t397) * sin(t218) + (pkin(1) * t113 + t131 * t256 + t259 * t490) * t512;
t222 = sin(t269);
t225 = cos(t269);
t280 = sin(qJ(1,2));
t286 = cos(qJ(1,2));
t125 = -t280 * t222 + t286 * t225;
t128 = t286 * t222 + t280 * t225;
t229 = 0.1e1 / t243;
t443 = t110 * t113;
t77 = (t125 * t291 + t128 * t290 + t289 * t443) * t229;
t563 = t110 / t315 * t56 * t77;
t281 = sin(qJ(2,1));
t418 = t265 * t281;
t378 = pkin(3) * t418;
t287 = cos(qJ(2,1));
t422 = t548 * t287;
t111 = 0.1e1 / (-t378 + t422);
t273 = -qJ(3,1) - pkin(5);
t244 = -pkin(6) + t273;
t316 = t244 ^ 2;
t114 = t281 * t548 + t287 * t486;
t159 = -t244 * t291 + t240;
t162 = t244 * t290 + t238;
t270 = legFrame(1,3);
t219 = t270 + qJ(1,1);
t299 = 0.2e1 * qJ(2,1);
t257 = sin(t299);
t260 = cos(t299);
t381 = 0.2e1 * t422;
t57 = (t260 * t384 + t162 * t381 + (-t162 * t418 - t257 * t426) * t393 + t396) * cos(t219) + (t260 * t385 + t159 * t381 + (-t159 * t418 - t257 * t427) * t393 + t397) * sin(t219) + (pkin(1) * t114 + t131 * t257 + t260 * t490) * t512;
t223 = sin(t270);
t226 = cos(t270);
t282 = sin(qJ(1,1));
t288 = cos(qJ(1,1));
t126 = -t282 * t223 + t288 * t226;
t129 = t288 * t223 + t282 * t226;
t231 = 0.1e1 / t244;
t441 = t111 * t114;
t78 = (t126 * t291 + t129 * t290 + t289 * t441) * t231;
t562 = t111 / t316 * t57 * t78;
t561 = 0.3e1 * t304 + 0.6e1 * t302;
t560 = 0.6e1 * t304 + 0.3e1 * t302;
t237 = t287 * pkin(2);
t254 = qJ(2,1) + pkin(7);
t216 = cos(t254);
t487 = pkin(3) * t216;
t555 = -t487 - t237;
t156 = 0.1e1 / t555 ^ 2;
t263 = t287 ^ 2;
t274 = xDDP(3);
t275 = xDDP(2);
t276 = xDDP(1);
t440 = t111 * t231;
t360 = t114 * t440;
t366 = t57 * t440;
t54 = t366 / 0.2e1;
t67 = t78 * pkin(1);
t41 = -t67 + t54;
t264 = t289 ^ 2;
t425 = t556 * t264;
t434 = t129 * t231;
t437 = t126 * t231;
t248 = t266 ^ 2;
t395 = t198 + t304 / 0.2e1;
t517 = -0.2e1 * (t248 - 0.1e1 / 0.2e1) * t302 - 0.2e1 * t395;
t531 = (t248 - 0.1e1) * t302;
t33 = -t276 * t437 - t275 * t434 - t274 * t360 - (t41 * t378 - (t263 * t517 + t531) * t78 + (-0.2e1 * t78 * t378 - t41) * t422) * t78 * t440 - t156 * t231 / (t237 + (t266 * t287 - t418) * pkin(3)) * t425 + t562 / 0.2e1;
t30 = t33 * pkin(1);
t239 = g(2) * t288;
t484 = g(1) * t282;
t172 = -t239 + t484;
t482 = g(2) * t282;
t483 = g(1) * t288;
t173 = t482 + t483;
t398 = -t172 * t226 - t173 * t223;
t428 = t156 * t264;
t559 = pkin(5) * t428 - 0.2e1 * t30 + t398;
t236 = t285 * pkin(2);
t252 = qJ(2,2) + pkin(7);
t214 = cos(t252);
t488 = pkin(3) * t214;
t554 = -t488 - t236;
t154 = 0.1e1 / t554 ^ 2;
t262 = t285 ^ 2;
t442 = t110 * t229;
t361 = t113 * t442;
t435 = t128 * t229;
t438 = t125 * t229;
t367 = t56 * t442;
t53 = t367 / 0.2e1;
t71 = pkin(1) * t77;
t45 = -t71 + t53;
t32 = -t276 * t438 - t275 * t435 - t274 * t361 - (t45 * t379 - (t262 * t517 + t531) * t77 + (-0.2e1 * t77 * t379 - t45) * t423) * t77 * t442 - t154 * t229 / (t236 + (t266 * t285 - t419) * pkin(3)) * t425 + t563 / 0.2e1;
t29 = t32 * pkin(1);
t234 = t280 * g(1);
t473 = t286 * g(2);
t170 = -t234 + t473;
t474 = t286 * g(1);
t477 = t280 * g(2);
t171 = t474 + t477;
t399 = t170 * t225 - t171 * t222;
t430 = t154 * t264;
t558 = pkin(5) * t430 - 0.2e1 * t29 + t399;
t235 = t283 * pkin(2);
t250 = qJ(2,3) + pkin(7);
t212 = cos(t250);
t489 = pkin(3) * t212;
t553 = -t489 - t235;
t152 = 0.1e1 / t553 ^ 2;
t261 = t283 ^ 2;
t444 = t109 * t227;
t362 = t112 * t444;
t368 = t55 * t444;
t52 = t368 / 0.2e1;
t69 = pkin(1) * t76;
t43 = -t69 + t52;
t436 = t127 * t227;
t439 = t124 * t227;
t31 = -t276 * t439 - t275 * t436 - t274 * t362 - (t43 * t380 - (t261 * t517 + t531) * t76 + (-0.2e1 * t76 * t380 - t43) * t424) * t76 * t444 - t152 * t227 / (t235 + (t266 * t283 - t420) * pkin(3)) * t425 + t564 / 0.2e1;
t28 = t31 * pkin(1);
t233 = t278 * g(1);
t475 = t284 * g(2);
t168 = -t233 + t475;
t476 = t284 * g(1);
t478 = t278 * g(2);
t169 = t476 + t478;
t400 = t168 * t224 - t169 * t221;
t432 = t152 * t264;
t557 = pkin(5) * t432 - 0.2e1 * t28 + t400;
t151 = 0.1e1 / t553;
t249 = t295 + pkin(7);
t211 = cos(t249);
t549 = 0.2e1 * pkin(1);
t402 = pkin(3) * t549;
t191 = 0.2e1 * t250;
t535 = t304 * t258 + t302 * cos(t191);
t341 = (t212 * t402 + t241 + (t283 * t549 + (t211 + t266) * t393) * pkin(2) + t535) * t564;
t205 = cos(t293 + qJ(2,3));
t208 = cos(-pkin(7) + qJ(2,3));
t509 = pkin(2) * pkin(3);
t388 = sin(t249) * t509;
t532 = t304 * t255 + t302 * sin(t191);
t323 = 0.2e1 * t388 + t532;
t301 = pkin(3) * t302;
t485 = pkin(3) * t304;
t348 = -0.2e1 * t301 - 0.4e1 * t485;
t389 = -0.2e1 * t485;
t491 = pkin(2) * t302;
t391 = -0.2e1 * t491;
t433 = t151 * t289;
t303 = pkin(2) * t304;
t514 = -0.2e1 * t303 - 0.4e1 * t491;
t515 = -0.4e1 * pkin(1) * (t499 + t395);
t371 = (t323 * t76 * t242 - (t205 * t391 + t208 * t389 + t212 * t348 + t283 * t514 + t515) * t433) * t152 * t289;
t495 = pkin(2) * t277;
t332 = pkin(3) * sin(t250) + t495;
t100 = t332 * t549 + t323;
t450 = t100 * t274;
t194 = t235 + pkin(1);
t118 = t194 * t278 + t242 * t284;
t122 = t194 * t284 - t242 * t278;
t82 = t118 * t224 + t122 * t221 + t127 * t489;
t454 = t82 * t275;
t79 = -t118 * t221 + t122 * t224 + t124 * t489;
t457 = t79 * t276;
t294 = 0.3e1 * qJ(2,3);
t305 = pkin(1) ^ 2;
t390 = -0.3e1 * t485;
t392 = -0.3e1 * t491;
t401 = -0.8e1 * t509;
t538 = -0.8e1 * t198 - 0.4e1 * t241;
t502 = (-(-0.4e1 * t388 - 0.2e1 * t532) * t242 * t433 + 0.4e1 * t553 * (pkin(1) * t52 - (t305 + t314) * t76) + (t301 * cos(0.3e1 * t250) + t489 * t560 + t303 * cos(t294) + t235 * t561 - (cos(t293 + t294) + t205) * t392 - (cos(t294 + pkin(7)) + t208) * t390) * t76 + (t211 * t401 - 0.4e1 * t535 + t538) * (-t69 + t368 / 0.4e1)) * t76;
t552 = (-(t450 / 0.2e1 + t502 / 0.4e1) * t151 + t457 + t454 - t371 / 0.2e1) * t227 + t28 + t151 * t341 / 0.4e1;
t153 = 0.1e1 / t554;
t251 = t297 + pkin(7);
t213 = cos(t251);
t192 = 0.2e1 * t252;
t536 = t304 * t259 + t302 * cos(t192);
t340 = (t214 * t402 + t241 + (t285 * t549 + (t213 + t266) * t393) * pkin(2) + t536) * t563;
t206 = cos(t293 + qJ(2,2));
t209 = cos(-pkin(7) + qJ(2,2));
t387 = sin(t251) * t509;
t533 = t304 * t256 + t302 * sin(t192);
t322 = 0.2e1 * t387 + t533;
t431 = t153 * t289;
t370 = (t322 * t77 * t243 - (t206 * t391 + t209 * t389 + t214 * t348 + t285 * t514 + t515) * t431) * t154 * t289;
t494 = pkin(2) * t279;
t331 = pkin(3) * sin(t252) + t494;
t101 = t331 * t549 + t322;
t448 = t101 * t274;
t195 = t236 + pkin(1);
t119 = t195 * t280 + t243 * t286;
t123 = t195 * t286 - t243 * t280;
t83 = t119 * t225 + t123 * t222 + t128 * t488;
t453 = t83 * t275;
t80 = -t119 * t222 + t123 * t225 + t125 * t488;
t456 = t80 * t276;
t296 = 0.3e1 * qJ(2,2);
t501 = (-(-0.4e1 * t387 - 0.2e1 * t533) * t243 * t431 + 0.4e1 * t554 * (pkin(1) * t53 - (t305 + t315) * t77) + (t301 * cos(0.3e1 * t252) + t488 * t560 + t303 * cos(t296) + t236 * t561 - (cos(t293 + t296) + t206) * t392 - (cos(t296 + pkin(7)) + t209) * t390) * t77 + (t213 * t401 - 0.4e1 * t536 + t538) * (-t71 + t367 / 0.4e1)) * t77;
t551 = (-t153 * (t448 / 0.2e1 + t501 / 0.4e1) + t456 + t453 - t370 / 0.2e1) * t229 + t29 + t153 * t340 / 0.4e1;
t155 = 0.1e1 / t555;
t253 = t299 + pkin(7);
t215 = cos(t253);
t193 = 0.2e1 * t254;
t537 = t304 * t260 + t302 * cos(t193);
t339 = (t216 * t402 + t241 + (t287 * t549 + (t215 + t266) * t393) * pkin(2) + t537) * t562;
t207 = cos(t293 + qJ(2,1));
t210 = cos(-pkin(7) + qJ(2,1));
t386 = sin(t253) * t509;
t534 = t304 * t257 + t302 * sin(t193);
t321 = 0.2e1 * t386 + t534;
t429 = t155 * t289;
t369 = (t321 * t78 * t244 - (t207 * t391 + t210 * t389 + t216 * t348 + t287 * t514 + t515) * t429) * t156 * t289;
t493 = pkin(2) * t281;
t330 = pkin(3) * sin(t254) + t493;
t102 = t330 * t549 + t321;
t446 = t102 * t274;
t196 = t237 + pkin(1);
t120 = t196 * t282 + t244 * t288;
t121 = t196 * t288 - t244 * t282;
t84 = t120 * t226 + t121 * t223 + t129 * t487;
t452 = t84 * t275;
t81 = -t120 * t223 + t121 * t226 + t126 * t487;
t455 = t81 * t276;
t298 = 0.3e1 * qJ(2,1);
t500 = (-(-0.4e1 * t386 - 0.2e1 * t534) * t244 * t429 + 0.4e1 * t555 * (pkin(1) * t54 - (t305 + t316) * t78) + (t301 * cos(0.3e1 * t254) + t487 * t560 + t303 * cos(t298) + t237 * t561 - (cos(t298 + t293) + t207) * t392 - (cos(t298 + pkin(7)) + t210) * t390) * t78 + (t215 * t401 - 0.4e1 * t537 + t538) * (-t67 + t366 / 0.4e1)) * t78;
t550 = (-t155 * (t446 / 0.2e1 + t500 / 0.4e1) + t455 + t452 - t369 / 0.2e1) * t231 + t30 + t155 * t339 / 0.4e1;
t541 = t168 * t221 + t169 * t224;
t540 = t170 * t222 + t171 * t225;
t539 = -t172 * t223 + t173 * t226;
t377 = t76 * t433;
t97 = (t332 * t432 + t274) * t151;
t508 = -t97 / 0.2e1;
t58 = pkin(1) * t377 + pkin(5) * t508;
t523 = -0.2e1 * t58;
t376 = t77 * t431;
t98 = (t331 * t430 + t274) * t153;
t507 = -t98 / 0.2e1;
t59 = pkin(1) * t376 + pkin(5) * t507;
t522 = -0.2e1 * t59;
t375 = t78 * t429;
t99 = (t330 * t428 + t274) * t155;
t506 = -t99 / 0.2e1;
t60 = pkin(1) * t375 + pkin(5) * t506;
t521 = -0.2e1 * t60;
t511 = -pkin(5) / 0.2e1;
t510 = pkin(1) * g(2);
t505 = pkin(5) * t31;
t504 = pkin(5) * t32;
t503 = pkin(5) * t33;
t498 = 0.2e1 * t261 - 0.1e1;
t497 = 0.2e1 * t262 - 0.1e1;
t496 = 0.2e1 * t263 - 0.1e1;
t492 = pkin(2) * t289;
t481 = g(3) * t283;
t480 = g(3) * t285;
t479 = g(3) * t287;
t466 = t271 * t31;
t465 = t272 * t32;
t464 = t273 * t33;
t463 = t283 * t31;
t462 = t285 * t32;
t461 = t287 * t33;
t460 = t31 * t277;
t459 = t32 * t279;
t458 = t33 * t281;
t451 = t100 * t151;
t449 = t101 * t153;
t447 = t102 * t155;
t417 = t277 * t283;
t416 = t279 * t285;
t415 = t281 * t287;
t73 = t76 ^ 2;
t374 = t73 * t417;
t74 = t77 ^ 2;
t373 = t74 * t416;
t75 = t78 ^ 2;
t372 = t75 * t415;
t365 = pkin(1) * t73 + t541;
t364 = pkin(1) * t74 + t540;
t363 = pkin(1) * t75 + t539;
t347 = 0.2e1 * t377;
t346 = 0.2e1 * t376;
t345 = 0.2e1 * t375;
t338 = t76 * t368;
t337 = t77 * t367;
t336 = t78 * t366;
t85 = -t277 * t97 + t283 * t432;
t87 = -t279 * t98 + t285 * t430;
t89 = -t281 * t99 + t287 * t428;
t317 = t151 * t460 + t153 * t459 + t155 * t458;
t300 = pkin(5) ^ 2;
t292 = pkin(1) * g(1);
t181 = -t273 * g(1) + t510;
t180 = t273 * g(2) + t292;
t179 = -t272 * g(1) + t510;
t178 = t272 * g(2) + t292;
t177 = -t271 * g(1) + t510;
t176 = t271 * g(2) + t292;
t149 = g(1) * t226 + g(2) * t223;
t148 = g(1) * t225 + g(2) * t222;
t147 = g(1) * t224 + g(2) * t221;
t146 = g(1) * t223 - g(2) * t226;
t145 = g(1) * t222 - g(2) * t225;
t144 = t221 * g(1) - t224 * g(2);
t108 = -t146 * t282 + t149 * t288;
t107 = -t145 * t280 + t148 * t286;
t106 = -t144 * t278 + t147 * t284;
t105 = t146 * t288 + t149 * t282;
t104 = t145 * t286 + t148 * t280;
t103 = t144 * t284 + t147 * t278;
t90 = -t281 * t428 - t287 * t99;
t88 = -t279 * t430 - t285 * t98;
t86 = -t277 * t432 - t283 * t97;
t27 = (t287 * t345 + t458) * t281;
t26 = (t285 * t346 + t459) * t279;
t25 = (t283 * t347 + t460) * t277;
t24 = t363 - t503;
t23 = t364 - t504;
t22 = t365 - t505;
t21 = t559 * t281 + t287 * t521;
t20 = t281 * t521 - t559 * t287;
t19 = t279 * t522 - t558 * t285;
t18 = t558 * t279 + t285 * t522;
t17 = t277 * t523 - t557 * t283;
t16 = t557 * t277 + t283 * t523;
t15 = t33 * t415 + t496 * t375;
t14 = t32 * t416 + t497 * t376;
t13 = t31 * t417 + t498 * t377;
t12 = -pkin(2) * t89 + t336 - 0.2e1 * t464 - t539;
t11 = -pkin(2) * t87 + t337 - 0.2e1 * t465 - t540;
t10 = -pkin(2) * t85 + t338 - 0.2e1 * t466 - t541;
t6 = t273 * t75 + (t281 * t345 - t461) * pkin(2) + t398 - t550;
t5 = t272 * t74 + (t279 * t346 - t462) * pkin(2) + t399 - t551;
t4 = t271 * t73 + (t277 * t347 - t463) * pkin(2) + t400 - t552;
t3 = 0.2e1 * ((-t239 / 0.2e1 + t484 / 0.2e1) * t226 + (t483 / 0.2e1 + t482 / 0.2e1) * t223 + t30 + (t511 - qJ(3,1) / 0.2e1) * t428 - (t281 * t78 * t492 - t339 / 0.8e1) * t155 - (-t455 / 0.2e1 - t452 / 0.2e1 + t369 / 0.4e1 - (-t446 / 0.4e1 - t500 / 0.8e1) * t155) * t231) * t237 + (t180 * t282 - t181 * t288) * t226 + (t180 * t288 + t181 * t282) * t223 - 0.2e1 * (qJ(3,1) * t506 + t60) * t493 + (t336 + 0.2e1 * t503) * qJ(3,1) + pkin(5) * t336 + pkin(1) * t550 + (qJ(3,1) ^ 2 + t263 * t304 + t300) * t33;
t2 = 0.2e1 * ((-t473 / 0.2e1 + t234 / 0.2e1) * t225 + (t474 / 0.2e1 + t477 / 0.2e1) * t222 + t29 + (t511 - qJ(3,2) / 0.2e1) * t430 - (t279 * t77 * t492 - t340 / 0.8e1) * t153 - (-t456 / 0.2e1 - t453 / 0.2e1 + t370 / 0.4e1 - (-t448 / 0.4e1 - t501 / 0.8e1) * t153) * t229) * t236 + (t178 * t280 - t179 * t286) * t225 + (t178 * t286 + t179 * t280) * t222 - 0.2e1 * (qJ(3,2) * t507 + t59) * t494 + (t337 + 0.2e1 * t504) * qJ(3,2) + pkin(5) * t337 + pkin(1) * t551 + (qJ(3,2) ^ 2 + t262 * t304 + t300) * t32;
t1 = 0.2e1 * ((-t475 / 0.2e1 + t233 / 0.2e1) * t224 + (t476 / 0.2e1 + t478 / 0.2e1) * t221 + t28 + (t511 - qJ(3,3) / 0.2e1) * t432 - (t277 * t76 * t492 - t341 / 0.8e1) * t151 - (-t457 / 0.2e1 - t454 / 0.2e1 + t371 / 0.4e1 - (-t450 / 0.4e1 - t502 / 0.8e1) * t151) * t227) * t235 + (t278 * t176 - t177 * t284) * t224 + (t176 * t284 + t278 * t177) * t221 - 0.2e1 * (qJ(3,3) * t508 + t58) * t495 + (t338 + 0.2e1 * t505) * qJ(3,3) + pkin(5) * t338 + pkin(1) * t552 + (qJ(3,3) ^ 2 + t261 * t304 + t300) * t31;
t7 = [-t31 * t439 - t32 * t438 - t33 * t437, -t103 * t439 - t104 * t438 - t105 * t437, -t106 * t439 - t107 * t438 - t108 * t437, -t25 * t439 - t26 * t438 - t27 * t437, -0.2e1 * t13 * t439 - 0.2e1 * t14 * t438 - 0.2e1 * t15 * t437, -t437 * t89 - t438 * t87 - t439 * t85, -t437 * t90 - t438 * t88 - t439 * t86, 0, -t17 * t439 - t19 * t438 - t20 * t437, -t16 * t439 - t18 * t438 - t21 * t437, -(t12 * t126 - t75 * t81) * t231 - (t11 * t125 - t74 * t80) * t229 - (t10 * t124 - t73 * t79) * t227, -(t126 * t3 + t6 * t81) * t231 - (t125 * t2 + t5 * t80) * t229 - (t1 * t124 + t4 * t79) * t227, t276 - g(1); -t31 * t436 - t32 * t435 - t33 * t434, -t103 * t436 - t104 * t435 - t105 * t434, -t106 * t436 - t107 * t435 - t108 * t434, -t25 * t436 - t26 * t435 - t27 * t434, -0.2e1 * t13 * t436 - 0.2e1 * t14 * t435 - 0.2e1 * t15 * t434, -t434 * t89 - t435 * t87 - t436 * t85, -t434 * t90 - t435 * t88 - t436 * t86, 0, -t17 * t436 - t19 * t435 - t20 * t434, -t16 * t436 - t18 * t435 - t21 * t434, -(t12 * t129 - t75 * t84) * t231 - (t11 * t128 - t74 * t83) * t229 - (t10 * t127 - t73 * t82) * t227, -(t129 * t3 + t6 * t84) * t231 - (t128 * t2 + t5 * t83) * t229 - (t1 * t127 + t4 * t82) * t227, t275 - g(2); -t31 * t362 - t32 * t361 - t33 * t360, -t103 * t362 - t104 * t361 - t105 * t360, -t106 * t362 - t107 * t361 - t108 * t360, t151 * t374 + t153 * t373 + t155 * t372 - t25 * t362 - t26 * t361 - t27 * t360, t151 * t498 * t73 + t153 * t497 * t74 + t155 * t496 * t75 - 0.2e1 * t13 * t362 - 0.2e1 * t14 * t361 - 0.2e1 * t15 * t360, -t360 * t89 - t361 * t87 - t362 * t85 - t317, -t151 * t463 - t153 * t462 - t155 * t461 - t360 * t90 - t361 * t88 - t362 * t86, t151 * t97 + t153 * t98 + t155 * t99, -t20 * t360 - t155 * (t24 * t281 - t479) - t19 * t361 - t153 * (t23 * t279 - t480) - t17 * t362 - t151 * (t22 * t277 - t481), -t21 * t360 - t155 * (g(3) * t281 + t24 * t287) - t18 * t361 - t153 * (g(3) * t279 + t23 * t285) - t16 * t362 - t151 * (g(3) * t277 + t22 * t283), -(t12 * t441 + t75 * t447 / 0.2e1) * t231 - (t11 * t443 + t74 * t449 / 0.2e1) * t229 - (t10 * t445 + t73 * t451 / 0.2e1) * t227 + t317 * pkin(2), -(t3 * t441 - t6 * t447 / 0.2e1) * t231 - (t2 * t443 - t5 * t449 / 0.2e1) * t229 - (t1 * t445 - t4 * t451 / 0.2e1) * t227 + (-t155 * ((t363 - t336 + t464) * t281 - t479) - t153 * ((t364 - t337 + t465) * t279 - t480) - t151 * ((t365 - t338 + t466) * t277 - t481) + (-t155 * (-t99 + t372) - t153 * (-t98 + t373) - t151 * (-t97 + t374)) * pkin(2)) * pkin(2), t274 - g(3);];
tauX_reg  = t7;
