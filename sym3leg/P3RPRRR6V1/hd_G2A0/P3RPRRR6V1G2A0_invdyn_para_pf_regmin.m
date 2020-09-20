% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR6V1G2A0
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
% tauX_reg [3x12]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR6V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:36:35
% EndTime: 2020-08-06 18:36:52
% DurationCPUTime: 16.28s
% Computational Cost: add. (116913->617), mult. (137418->1031), div. (18858->22), fcn. (98790->102), ass. (0->446)
t261 = sin(pkin(7));
t289 = -pkin(6) - pkin(5);
t170 = t289 * t261;
t157 = -t170 + pkin(1);
t273 = sin(qJ(1,3));
t126 = t157 * t273;
t262 = cos(pkin(7));
t279 = cos(qJ(1,3));
t445 = t279 * t289;
t452 = t261 * t279;
t579 = t126 + pkin(2) * t452 + (pkin(2) * t273 + t445) * t262;
t275 = sin(qJ(1,2));
t127 = t157 * t275;
t281 = cos(qJ(1,2));
t444 = t281 * t289;
t451 = t261 * t281;
t578 = t127 + pkin(2) * t451 + (pkin(2) * t275 + t444) * t262;
t277 = sin(qJ(1,1));
t128 = t157 * t277;
t283 = cos(qJ(1,1));
t443 = t283 * t289;
t450 = t261 * t283;
t577 = t128 + pkin(2) * t450 + (pkin(2) * t277 + t443) * t262;
t311 = pkin(2) ^ 2;
t291 = 0.3e1 * t311;
t309 = pkin(3) ^ 2;
t576 = t291 + 0.3e1 / 0.4e1 * t309;
t266 = legFrame(3,2);
t202 = sin(t266);
t205 = cos(t266);
t123 = g(1) * t205 - t202 * g(2);
t100 = g(3) * t279 + t123 * t273;
t229 = qJ(1,3) + pkin(7);
t161 = t266 + t229;
t162 = -t266 + t229;
t107 = -sin(t161) + sin(t162);
t110 = cos(t162) + cos(t161);
t278 = cos(qJ(3,3));
t208 = t278 * pkin(3);
t177 = t208 + pkin(2);
t199 = t262 * pkin(1);
t141 = t199 + t177;
t129 = 0.1e1 / t141;
t130 = 0.1e1 / t141 ^ 2;
t189 = sin(t229);
t270 = xDDP(2);
t271 = xDDP(1);
t272 = sin(qJ(3,3));
t502 = pkin(3) * t272;
t169 = t289 * t262;
t286 = xDP(3);
t287 = xDP(2);
t288 = xDP(1);
t310 = 0.1e1 / pkin(3);
t398 = pkin(2) * t262 + t157;
t436 = pkin(2) * t261 + t169;
t167 = t199 + pkin(2);
t448 = t272 * t167;
t254 = 0.1e1 / t272;
t468 = t129 * t254;
t300 = 0.2e1 * qJ(3,3);
t238 = sin(t300);
t505 = pkin(3) * t238;
t53 = (-t286 * ((t262 * t208 + t398) * t279 - t273 * (t261 * t208 + t436)) / (t505 / 0.2e1 + t448) + (t202 * t287 - t205 * t288) * ((t177 * t262 + t157) * t273 + t279 * (t177 * t261 + t169)) * t468) * t310;
t415 = t53 * t502;
t509 = pkin(1) * t261;
t336 = -t289 / 0.2e1 + t509 / 0.2e1;
t521 = -0.2e1 * t286;
t74 = t107 * t287 + t110 * t288 + t189 * t521;
t43 = t129 * t336 * t74 - t415;
t269 = xDDP(3);
t467 = t129 * t269;
t491 = t74 / 0.2e1;
t40 = -t189 * t467 + (-t43 + t415) * t130 * t491 + (t107 * t270 + t110 * t271) * t129 / 0.2e1;
t50 = t53 ^ 2;
t516 = pkin(1) * t40;
t544 = g(3) * t273 - t123 * t279;
t561 = 0.2e1 * pkin(2);
t575 = (t544 + 0.2e1 * t516) * t262 - t50 * pkin(5) - (pkin(1) * t50 - t100) * t261 + t40 * t561;
t267 = legFrame(2,2);
t203 = sin(t267);
t206 = cos(t267);
t124 = g(1) * t206 - t203 * g(2);
t101 = g(3) * t281 + t124 * t275;
t232 = qJ(1,2) + pkin(7);
t163 = t267 + t232;
t164 = -t267 + t232;
t108 = -sin(t163) + sin(t164);
t111 = cos(t164) + cos(t163);
t280 = cos(qJ(3,2));
t209 = t280 * pkin(3);
t178 = t209 + pkin(2);
t142 = t199 + t178;
t132 = 0.1e1 / t142;
t133 = 0.1e1 / t142 ^ 2;
t190 = sin(t232);
t274 = sin(qJ(3,2));
t501 = pkin(3) * t274;
t447 = t274 * t167;
t303 = 0.2e1 * qJ(3,2);
t241 = sin(t303);
t504 = pkin(3) * t241;
t377 = t286 * ((t262 * t209 + t398) * t281 - t275 * (t261 * t209 + t436)) / (t504 / 0.2e1 + t447);
t255 = 0.1e1 / t274;
t465 = t132 * t255;
t391 = ((t178 * t262 + t157) * t275 + t281 * (t178 * t261 + t169)) * t465;
t456 = t206 * t288;
t457 = t203 * t287;
t54 = (-t377 + (-t456 + t457) * t391) * t310;
t414 = t54 * t501;
t75 = t108 * t287 + t111 * t288 + t190 * t521;
t551 = t132 * t75;
t44 = t336 * t551 - t414;
t464 = t132 * t269;
t490 = t75 / 0.2e1;
t41 = -t190 * t464 + (-t44 + t414) * t133 * t490 + (t108 * t270 + t111 * t271) * t132 / 0.2e1;
t51 = t54 ^ 2;
t515 = pkin(1) * t41;
t543 = g(3) * t275 - t124 * t281;
t574 = (t543 + 0.2e1 * t515) * t262 - t51 * pkin(5) - (pkin(1) * t51 - t101) * t261 + t41 * t561;
t268 = legFrame(1,2);
t204 = sin(t268);
t207 = cos(t268);
t125 = g(1) * t207 - t204 * g(2);
t102 = g(3) * t283 + t125 * t277;
t235 = qJ(1,1) + pkin(7);
t165 = t268 + t235;
t166 = -t268 + t235;
t109 = -sin(t165) + sin(t166);
t112 = cos(t166) + cos(t165);
t282 = cos(qJ(3,1));
t210 = t282 * pkin(3);
t180 = t210 + pkin(2);
t143 = t199 + t180;
t135 = 0.1e1 / t143;
t136 = 0.1e1 / t143 ^ 2;
t191 = sin(t235);
t276 = sin(qJ(3,1));
t500 = pkin(3) * t276;
t256 = 0.1e1 / t276;
t462 = t135 * t256;
t373 = t207 * t462;
t446 = t276 * t167;
t306 = 0.2e1 * qJ(3,1);
t244 = sin(t306);
t503 = pkin(3) * t244;
t87 = (t180 * t262 + t157) * t277 + t283 * (t180 * t261 + t169);
t545 = (0.1e1 / (t503 / 0.2e1 + t446) * t286 * ((t262 * t210 + t398) * t283 - t277 * (t261 * t210 + t436)) + t288 * t87 * t373) * t310;
t371 = t310 * t462;
t70 = t204 * t287 * t87 * t371;
t55 = t70 - t545;
t413 = t55 * t500;
t76 = t109 * t287 + t112 * t288 + t191 * t521;
t549 = t135 * t76;
t45 = t336 * t549 - t413;
t461 = t135 * t269;
t489 = t76 / 0.2e1;
t42 = -t191 * t461 + (-t45 + t413) * t136 * t489 + (t109 * t270 + t112 * t271) * t135 / 0.2e1;
t514 = pkin(1) * t42;
t52 = t55 ^ 2;
t542 = g(3) * t277 - t125 * t283;
t573 = (t542 + 0.2e1 * t514) * t262 - t52 * pkin(5) - (pkin(1) * t52 - t102) * t261 + t42 * t561;
t423 = 0.4e1 * pkin(3);
t572 = -0.2e1 * t545 + 0.2e1 * t70;
t571 = (t377 / 0.6e1 + (t456 / 0.6e1 - t457 / 0.6e1) * t391) * t310;
t570 = (t377 / 0.4e1 + (t456 / 0.4e1 - t457 / 0.4e1) * t391) * t309 * t310;
t569 = (t377 / 0.3e1 + (t456 / 0.3e1 - t457 / 0.3e1) * t391) * t310;
t568 = (t377 / 0.2e1 + (t456 / 0.2e1 - t457 / 0.2e1) * t391) * t310;
t39 = t542 + t514;
t38 = t543 + t515;
t37 = t544 + t516;
t560 = 0.4e1 * pkin(2);
t424 = 0.2e1 * pkin(3);
t520 = 0.2e1 * t289;
t472 = t112 * t135;
t558 = t472 / 0.2e1;
t474 = t111 * t132;
t557 = t474 / 0.2e1;
t476 = t110 * t129;
t556 = t476 / 0.2e1;
t478 = t109 * t135;
t555 = t478 / 0.2e1;
t480 = t108 * t132;
t554 = t480 / 0.2e1;
t482 = t107 * t129;
t553 = t482 / 0.2e1;
t460 = t177 * t273;
t79 = (t445 + t460) * t262 + t126 + t177 * t452;
t552 = t129 * t79;
t459 = t178 * t275;
t81 = (t444 + t459) * t262 + t127 + t178 * t451;
t550 = t132 * t81;
t458 = t180 * t277;
t83 = (t443 + t458) * t262 + t128 + t180 * t450;
t494 = t135 * t83;
t548 = t130 * t278;
t547 = t133 * t280;
t546 = t136 * t282;
t297 = 0.2e1 * pkin(7);
t218 = cos(t297);
t312 = pkin(1) ^ 2;
t184 = t312 * t218;
t538 = -t184 - t309 / 0.2e1;
t532 = -0.2e1 * pkin(1);
t531 = -0.2e1 * pkin(2);
t416 = pkin(2) * t199;
t292 = 0.2e1 * t311;
t429 = t292 + t312;
t119 = 0.4e1 * t416 + t184 + t429;
t221 = pkin(7) + qJ(3,2);
t171 = 0.2e1 * t221;
t226 = -pkin(7) + qJ(3,2);
t173 = 0.2e1 * t226;
t185 = sin(t221);
t187 = sin(t226);
t220 = pkin(7) + t303;
t192 = cos(t220);
t193 = cos(t221);
t225 = -pkin(7) + t303;
t195 = cos(t225);
t196 = cos(t226);
t302 = 0.3e1 * qJ(3,2);
t219 = pkin(7) + t302;
t224 = -pkin(7) + t302;
t230 = qJ(3,2) + t297;
t231 = qJ(3,2) - 0.2e1 * pkin(7);
t301 = 0.4e1 * qJ(3,2);
t239 = sin(t301);
t240 = sin(t302);
t249 = cos(t302);
t250 = cos(t303);
t296 = -0.2e1 * t312;
t308 = pkin(3) * t309;
t407 = pkin(1) * t169;
t217 = sin(t297);
t437 = t312 * t217;
t410 = 0.2e1 * t437;
t332 = -0.4e1 * t407 + t410;
t260 = t289 ^ 2;
t68 = t132 * t490;
t67 = t312 * t68;
t343 = t67 + (t260 + t576) * t68;
t435 = t260 / 0.2e1 + t312;
t345 = 0.3e1 / 0.8e1 * t309 + t311 / 0.2e1 + t435;
t346 = -0.2e1 * t184 - 0.4e1 * t311 + t296 - t309;
t349 = t274 * t68;
t160 = -t289 + t509;
t363 = t160 * t424;
t368 = -t310 / 0.4e1;
t369 = -0.2e1 * t311 + t538;
t392 = t81 * t465;
t513 = pkin(3) * t54;
t428 = 0.4e1 * t513;
t431 = t249 - t280;
t440 = t309 * cos(t301);
t449 = t269 * t310;
t46 = t68 + t568;
t47 = t68 - t568;
t498 = (t190 * t520 + cos(t232) * t531 + t281 * t532 + (-cos(qJ(1,2) - t226) - cos(qJ(1,2) + t221)) * pkin(3)) / (t274 * t561 + t504 + (t185 + t187) * pkin(1));
t506 = pkin(2) * t160;
t507 = pkin(1) * t289;
t517 = pkin(2) * pkin(3);
t519 = -0.6e1 * t309;
t525 = -0.2e1 * t262;
t526 = -0.2e1 * t309 - 0.2e1 * t119;
t361 = -0.2e1 * t407 + t437;
t96 = t361 + 0.2e1 * t506;
t12 = t449 * t498 + (t309 * t54 * t249 * t520 + t410 * t513 + (t506 + (-t199 - t209 / 0.2e1) * t289) * t428 + (-0.8e1 * (t309 / 0.4e1 + 0.3e1 / 0.2e1 * t311 + t435) * t504 - t308 * t239) * t68 + (t46 * sin(t173) + t47 * sin(t171)) * pkin(3) * t296 + (-0.4e1 * (t343 - t570) * t187 - 0.4e1 * (t343 + t570) * t185 - 0.3e1 * ((t68 + t569) * sin(t224) + (t68 - t569) * sin(t219)) * t309 - 0.12e2 * ((t68 + t571) * sin(t225) + (t68 - t571) * sin(t220)) * t517) * pkin(1) + (-t46 * t192 + t47 * t195) * t507 * t423 + 0.4e1 * (cos(t231) - cos(t230)) * t289 * t67 + ((t240 * t519 + 0.8e1 * (-t193 + t196) * t507) * t68 - 0.4e1 * (sin(t231) + sin(t230)) * t67 + t289 * t250 * t428 - 0.16e2 * t345 * t349) * pkin(2)) / (t292 * t250 + t440 / 0.2e1 + (cos(t173) / 0.2e1 + cos(t171) / 0.2e1 + t250 - 0.1e1) * t312 + t431 * pkin(2) * t424 + ((t192 + t195 + t525) * t561 + (cos(t224) + cos(t219) - t193 - t196) * pkin(3)) * pkin(1) + t369) * t368 * t551 - (t44 * t560 + (t241 * t526 - t309 * t239 + (-t167 * t240 - t274 * t199) * t423) * t54 + (-0.2e1 * t250 * t96 - t431 * t363 + t332) * t68) / (0.2e1 * t119 * t250 + t440 - 0.4e1 * (t209 + t561) * t199 + (-pkin(2) * t280 + t167 * t249) * t423 + t346) * t54 + (t203 * t270 - t206 * t271) * t310 * t392;
t358 = t54 * t68;
t510 = t261 / 0.2e1;
t518 = pkin(5) / 0.2e1;
t529 = -0.2e1 * pkin(2) * t358 - 0.2e1 * t12 * t518 - 0.2e1 * (t12 * t510 + t262 * t358) * pkin(1);
t223 = pkin(7) + qJ(3,1);
t172 = 0.2e1 * t223;
t228 = -pkin(7) + qJ(3,1);
t174 = 0.2e1 * t228;
t179 = t210 + t561;
t186 = sin(t223);
t188 = sin(t228);
t222 = pkin(7) + t306;
t227 = -pkin(7) + t306;
t305 = 0.3e1 * qJ(3,1);
t233 = t305 + pkin(7);
t234 = t305 - pkin(7);
t304 = 0.4e1 * qJ(3,1);
t242 = sin(t304);
t243 = sin(t305);
t252 = cos(t305);
t253 = cos(t306);
t295 = 0.2e1 * t312;
t408 = pkin(1) * t170;
t430 = t260 + t312;
t322 = -0.8e1 * (t430 + t576) * t199 - 0.8e1 * (t345 - t408) * t561 + 0.8e1 * (-pkin(2) * t218 + t217 * t289) * t312;
t344 = -t312 + t369;
t364 = pkin(3) * (t252 - t282);
t370 = cos(t222) + cos(t227) + t525;
t411 = t167 * t519;
t412 = -0.2e1 * t160 * t309;
t158 = -0.2e1 * t408;
t427 = -0.4e1 * pkin(3) * (0.6e1 * t416 + t158 + t295 + t291 + t260 - t538);
t439 = t309 * cos(t304);
t69 = t135 * t489;
t48 = t69 - t572;
t49 = t69 + t572;
t497 = (t191 * t520 + cos(t235) * t531 + t283 * t532 + (-cos(qJ(1,1) - t228) - cos(qJ(1,1) + t223)) * pkin(3)) / (t276 * t561 + t503 + (t186 + t188) * pkin(1));
t15 = t449 * t497 + ((t252 * t412 + (t160 * t179 - t253 * t96 + t361) * t424) * t55 + (-t242 * t308 + t243 * t411 + t244 * t427 + t276 * t322) * t69) / (t119 * t253 + t439 / 0.2e1 - 0.2e1 * t179 * t199 + (-pkin(2) * t282 + t167 * t252) * t424 + t344) * t368 * t549 - (t45 * t560 + ((t69 - t55) * sin(t174) - (t69 + t55) * sin(t172)) * t312 + (-0.2e1 * (t309 + t429) * t244 - 0.4e1 * t243 * t517 - t309 * t242) * t55 + (t410 + (t253 * t560 + 0.2e1 * t364) * t289) * t69 + ((sin(t227) * t48 - sin(t222) * t49) * t561 + t370 * t69 * t520 + ((-sin(t233) - t188) * t49 + (sin(t234) + t186) * t48) * pkin(3)) * pkin(1)) / ((t295 + 0.4e1 * t311) * t253 + t439 + (cos(t174) + cos(t172)) * t312 + t364 * t560 + (t370 * t560 + (cos(t234) + cos(t233) - cos(t228) - cos(t223)) * t424) * pkin(1) + t346) * t55 + (t204 * t270 - t207 * t271) * t83 * t371;
t356 = t55 * t69;
t528 = -0.2e1 * pkin(2) * t356 - 0.2e1 * t15 * t518 - 0.2e1 * (t15 * t510 + t262 * t356) * pkin(1);
t176 = t208 + t561;
t298 = 0.4e1 * qJ(3,3);
t236 = sin(t298);
t299 = 0.3e1 * qJ(3,3);
t237 = sin(t299);
t246 = cos(t299);
t328 = pkin(3) * (-pkin(2) * t278 + t167 * t246);
t395 = t129 * t491;
t396 = t79 * t468;
t409 = t176 * t199;
t441 = t309 * cos(t298);
t247 = cos(t300);
t470 = t119 * t247;
t492 = t247 * t96;
t418 = pkin(7) + qJ(3,3);
t419 = -pkin(7) + qJ(3,3);
t499 = (t189 * t520 + cos(t229) * t531 + t279 * t532 + (-cos(qJ(1,3) - t419) - cos(qJ(1,3) + t418)) * pkin(3)) / (t272 * t561 + t505 + (sin(t418) + sin(t419)) * pkin(1));
t18 = -(t43 * t560 + (t238 * t526 - t309 * t236 + (-t167 * t237 - t272 * t199) * t423) * t53 + (-0.2e1 * t492 + (-t246 + t278) * t363 + t332) * t395) / (0.4e1 * t328 + t346 - 0.4e1 * t409 + t441 + 0.2e1 * t470) * t53 + (t269 * t499 - ((t246 * t412 + (t160 * t176 + t361 - t492) * t424) * t53 + (-t236 * t308 + t237 * t411 + t238 * t427 + t272 * t322) * t395) / (t470 + t441 / 0.2e1 - 0.2e1 * t409 + 0.2e1 * t328 + t344) * t395 / 0.2e1 + (t202 * t270 - t205 * t271) * t396) * t310;
t360 = t53 * t395;
t527 = -0.2e1 * pkin(2) * t360 - 0.2e1 * t18 * t518 - 0.2e1 * (t18 * t510 + t262 * t360) * pkin(1);
t524 = 0.2e1 * t278;
t523 = 0.2e1 * t280;
t522 = 0.2e1 * t282;
t512 = t53 * pkin(3);
t511 = t55 * pkin(3);
t508 = pkin(1) / 0.2e1;
t71 = t74 ^ 2;
t496 = t130 * t71;
t72 = t75 ^ 2;
t495 = t133 * t72;
t73 = t76 ^ 2;
t493 = t136 * t73;
t488 = t310 * t79;
t487 = t310 * t81;
t486 = t310 * t83;
t485 = t40 * t272;
t484 = t41 * t274;
t483 = t42 * t276;
t481 = t107 / 0.2e1;
t479 = t108 / 0.2e1;
t477 = t109 / 0.2e1;
t475 = t110 / 0.2e1;
t473 = t111 / 0.2e1;
t471 = t112 / 0.2e1;
t469 = t129 * t189;
t466 = t132 * t190;
t463 = t135 * t191;
t455 = t254 * t278;
t454 = t255 * t280;
t453 = t256 * t282;
t442 = t310 / 0.4e1;
t421 = t167 * t424;
t257 = t278 ^ 2;
t406 = pkin(3) * (t262 * t273 + t452) * t257;
t258 = t280 ^ 2;
t405 = pkin(3) * (t262 * t275 + t451) * t258;
t259 = t282 ^ 2;
t404 = pkin(3) * (t262 * t277 + t450) * t259;
t403 = t40 * t552;
t402 = t41 * t550;
t401 = t278 * t499;
t400 = t280 * t498;
t399 = t282 * t497;
t103 = t158 + t311 + 0.2e1 * t416 + t430;
t59 = -t202 * t406 + (-t579 * t202 + t205 * t502) * t278 + t205 * t448;
t60 = t205 * t406 + (t202 * t502 + t579 * t205) * t278 + t202 * t448;
t80 = (t177 * t279 - t273 * t289) * t262 + t157 * t279 - t261 * t460;
t28 = -t202 * g(1) - t205 * g(2) + (t80 * t278 * t467 + (-t160 * t415 + (t257 * t309 + t278 * t421 + t103) * t395) * t491 * t548) * t254 + (t60 * t271 + t59 * t270 - ((t160 * t272 * t395 - t512) * t278 - t53 * t167) * t512) * t468;
t397 = t28 * t468;
t394 = t496 / 0.4e1;
t61 = -t203 * t405 + (-t578 * t203 + t206 * t501) * t280 + t206 * t447;
t62 = t206 * t405 + (t203 * t501 + t578 * t206) * t280 + t203 * t447;
t82 = (t178 * t281 - t275 * t289) * t262 + t157 * t281 - t261 * t459;
t29 = -t203 * g(1) - t206 * g(2) + (t82 * t280 * t464 + (-t160 * t414 + (t258 * t309 + t280 * t421 + t103) * t68) * t490 * t547) * t255 + (t62 * t271 + t61 * t270 - ((t160 * t349 - t513) * t280 - t54 * t167) * t513) * t465;
t393 = t29 * t465;
t390 = t495 / 0.4e1;
t389 = t204 * t494;
t63 = -t204 * t404 + (-t577 * t204 + t207 * t500) * t282 + t207 * t446;
t64 = t207 * t404 + (t204 * t500 + t577 * t207) * t282 + t204 * t446;
t84 = (t180 * t283 - t277 * t289) * t262 + t157 * t283 - t261 * t458;
t30 = -t204 * g(1) - t207 * g(2) + (t84 * t282 * t461 + (-t160 * t413 + (t259 * t309 + t282 * t421 + t103) * t69) * t489 * t546) * t256 + (t64 * t271 + t63 * t270 - ((t160 * t276 * t69 - t511) * t282 - t55 * t167) * t511) * t462;
t388 = t30 * t462;
t387 = t493 / 0.4e1;
t386 = t202 * t488;
t385 = t203 * t487;
t384 = t204 * t486;
t383 = t205 * t488;
t382 = t206 * t487;
t381 = t207 * t486;
t380 = t80 * t455;
t379 = t82 * t454;
t378 = t84 * t453;
t375 = t129 * t455;
t374 = t132 * t454;
t372 = t135 * t453;
t359 = t71 * t548 * t552;
t357 = t72 * t547 * t550;
t355 = t73 * t546 * t494;
t354 = t202 * t396;
t353 = t205 * t396;
t352 = t203 * t392;
t351 = t206 * t392;
t348 = t256 * t389;
t347 = t83 * t373;
t339 = t40 * t79 * t375;
t338 = t41 * t81 * t374;
t337 = t42 * t83 * t372;
t58 = (-0.2e1 * t259 + 0.1e1) * t387;
t57 = (-0.2e1 * t258 + 0.1e1) * t390;
t56 = (-0.2e1 * t257 + 0.1e1) * t394;
t36 = (t356 * t522 + t483) * t276;
t35 = (t358 * t523 + t484) * t274;
t34 = (t360 * t524 + t485) * t272;
t33 = t483 * t522 + (0.4e1 * t259 - 0.2e1) * t356;
t32 = t484 * t523 + (0.4e1 * t258 - 0.2e1) * t358;
t31 = t485 * t524 + (0.4e1 * t257 - 0.2e1) * t360;
t27 = (pkin(1) * t387 + t102) * t262 - t39 * t261 + pkin(2) * t387 - t42 * pkin(5);
t26 = (pkin(1) * t390 + t101) * t262 - t38 * t261 + pkin(2) * t390 - t41 * pkin(5);
t25 = (pkin(1) * t394 + t100) * t262 - t37 * t261 + pkin(2) * t394 - t40 * pkin(5);
t24 = t27 * t282 - t276 * t30;
t23 = t27 * t276 + t282 * t30;
t22 = t26 * t280 - t274 * t29;
t21 = t26 * t274 + t280 * t29;
t20 = t25 * t278 - t272 * t28;
t19 = t25 * t272 + t278 * t28;
t17 = t18 * t278 - t272 * t50;
t16 = t18 * t272 + t278 * t50;
t14 = t15 * t282 - t276 * t52;
t13 = t15 * t276 + t282 * t52;
t11 = t12 * t280 - t274 * t51;
t10 = t12 * t274 + t280 * t51;
t8 = t272 * t527 + t575 * t278;
t7 = -t575 * t272 + t278 * t527;
t5 = t276 * t528 + t573 * t282;
t4 = -t573 * t276 + t282 * t528;
t2 = t274 * t529 + t574 * t280;
t1 = -t574 * t274 + t280 * t529;
t3 = [t40 * t556 + t41 * t557 + t42 * t558, t542 * t558 + t543 * t557 + t544 * t556, t100 * t556 + t101 * t557 + t102 * t558, t60 * t397 + t62 * t393 + t64 * t388 + (t37 * t476 + t38 * t474 + t39 * t472) * t508, (t205 * t359 + t206 * t357 + t207 * t355) * t442 + t34 * t556 + t35 * t557 + t36 * t558, (-t347 * t58 - t351 * t57 - t353 * t56) * t310 + t31 * t556 + t32 * t557 + t33 * t558, (-t207 * t42 * t494 - t205 * t403 - t206 * t402) * t310 + t10 * t557 + t13 * t558 + t16 * t556, (-t205 * t339 - t206 * t338 - t207 * t337) * t310 + t11 * t557 + t14 * t558 + t17 * t556, (-t12 * t351 - t15 * t347 - t18 * t353) * t310, (t5 * t471 + (t14 * t64 - t23 * t381) * t256) * t135 + (t2 * t473 + (t11 * t62 - t21 * t382) * t255) * t132 + (t8 * t475 + (t17 * t60 - t19 * t383) * t254) * t129, (t4 * t471 + (-t13 * t64 - t24 * t381) * t256) * t135 + (t1 * t473 + (-t10 * t62 - t22 * t382) * t255) * t132 + (t7 * t475 + (-t16 * t60 - t20 * t383) * t254) * t129, t271 - g(1); t40 * t553 + t41 * t554 + t42 * t555, t542 * t555 + t543 * t554 + t544 * t553, t100 * t553 + t101 * t554 + t102 * t555, t59 * t397 + t61 * t393 + t63 * t388 + (t37 * t482 + t38 * t480 + t39 * t478) * t508, (-t202 * t359 - t203 * t357 - t204 * t355) * t442 + t34 * t553 + t35 * t554 + t36 * t555, (t348 * t58 + t352 * t57 + t354 * t56) * t310 + t31 * t553 + t32 * t554 + t33 * t555, (t202 * t403 + t203 * t402 + t389 * t42) * t310 + t10 * t554 + t13 * t555 + t16 * t553, (t202 * t339 + t203 * t338 + t204 * t337) * t310 + t11 * t554 + t14 * t555 + t17 * t553, (t12 * t352 + t15 * t348 + t18 * t354) * t310, (t5 * t477 + (t14 * t63 + t23 * t384) * t256) * t135 + (t2 * t479 + (t11 * t61 + t21 * t385) * t255) * t132 + (t8 * t481 + (t17 * t59 + t19 * t386) * t254) * t129, (t4 * t477 + (-t13 * t63 + t24 * t384) * t256) * t135 + (t1 * t479 + (-t10 * t61 + t22 * t385) * t255) * t132 + (t7 * t481 + (-t16 * t59 + t20 * t386) * t254) * t129, t270 - g(2); -t40 * t469 - t41 * t466 - t42 * t463, -t463 * t542 - t466 * t543 - t469 * t544, -t100 * t469 - t101 * t466 - t102 * t463, t28 * t80 * t375 + t29 * t82 * t374 + t30 * t84 * t372 + (-t37 * t469 - t38 * t466 - t39 * t463) * pkin(1), -t34 * t469 - t35 * t466 - t36 * t463 + (-t272 * t401 * t496 - t274 * t400 * t495 - t276 * t399 * t493) * t442, -t31 * t469 - t32 * t466 - t33 * t463 + (t58 * t497 + t57 * t498 + t56 * t499) * t310, -t10 * t466 - t16 * t469 - t13 * t463 + (t483 * t497 + t484 * t498 + t485 * t499) * t310, -t11 * t466 - t17 * t469 - t14 * t463 + (t399 * t42 + t40 * t401 + t400 * t41) * t310, (t12 * t498 + t15 * t497 + t18 * t499) * t310, (t14 * t378 - t191 * t5) * t135 + (t11 * t379 - t190 * t2) * t132 + (t17 * t380 - t189 * t8) * t129 + (t19 * t499 + t21 * t498 + t23 * t497) * t310, (-t13 * t378 - t191 * t4) * t135 + (-t1 * t190 - t10 * t379) * t132 + (-t16 * t380 - t189 * t7) * t129 + (t20 * t499 + t22 * t498 + t24 * t497) * t310, t269 - g(3);];
tauX_reg  = t3;
