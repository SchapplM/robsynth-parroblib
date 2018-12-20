% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRP2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x11]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRP2A0_invdyn_para_pf_reg(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_reg: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_reg: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_reg: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_reg: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_reg: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_reg: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_reg: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_reg: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:39:07
% EndTime: 2018-12-20 17:39:15
% DurationCPUTime: 9.14s
% Computational Cost: add. (78519->515), mult. (156556->772), div. (2844->6), fcn. (50847->14), ass. (0->366)
t266 = (pkin(2) ^ 2);
t422 = -t266 - 1;
t240 = legFrame(3,3);
t218 = sin(t240);
t221 = cos(t240);
t170 = g(1) * t218 - g(2) * t221;
t246 = sin(qJ(2,3));
t249 = cos(qJ(2,3));
t256 = (qJ(3,3) ^ 2);
t212 = -t256 - t422;
t235 = t249 ^ 2;
t339 = t246 * t249;
t315 = pkin(2) * t339;
t145 = 0.2e1 * qJ(3,3) * t315 + t212 * t235 - t256 + t422;
t136 = 0.1e1 / t145;
t367 = qJ(3,3) * t218;
t194 = pkin(2) * t367;
t150 = t212 * t221 - 0.2e1 * t194;
t366 = qJ(3,3) * t221;
t310 = pkin(2) * t366;
t153 = t212 * t218 + 0.2e1 * t310;
t106 = t150 * t235 + t153 * t339 - t221 * t266 + t194 - t221;
t109 = -t150 * t339 + t153 * t235 - t218 * t266 - t218 - t310;
t255 = xP(3);
t227 = sin(t255);
t228 = cos(t255);
t259 = koppelP(3,2);
t262 = koppelP(3,1);
t164 = t227 * t262 + t228 * t259;
t167 = -t227 * t259 + t228 * t262;
t252 = xDP(3);
t238 = t252 ^ 2;
t243 = xDDP(3);
t244 = xDDP(2);
t124 = -t164 * t238 + t167 * t243 + t244;
t245 = xDDP(1);
t127 = -t164 * t243 - t167 * t238 + t245;
t300 = t106 * t127 + t109 * t124;
t265 = pkin(2) * t266;
t137 = 0.1e1 / t145 ^ 2;
t362 = qJ(3,3) * t259;
t188 = pkin(2) * t262 + t362;
t361 = qJ(3,3) * t262;
t189 = pkin(2) * t259 - t361;
t112 = (t188 * t228 - t189 * t227) * t221 + (t188 * t227 + t189 * t228) * t218;
t253 = xDP(2);
t226 = pkin(2) * t253;
t254 = xDP(1);
t363 = qJ(3,3) * t254;
t364 = qJ(3,3) * t253;
t421 = pkin(2) * t254;
t103 = (t226 - t363) * t221 + t218 * (-t364 - t421) + t112 * t252;
t291 = t164 * t218 + t167 * t221;
t365 = qJ(3,3) * t249;
t325 = -0.2e1 * t365;
t94 = t103 * t246 + (-t218 * t254 + t253 * t221 + t291 * t252) * t325;
t399 = t137 * t94;
t229 = t256 * t259;
t316 = pkin(2) * t361;
t173 = t229 + t259 + t316;
t215 = pkin(2) * t362;
t336 = t256 * t262;
t174 = t215 - t262 - t336;
t276 = (t173 * t227 + t174 * t228) * t221 - t218 * (t173 * t228 - t174 * t227);
t79 = ((-pkin(2) * t363 - t253 * t256 - t253) * t221 - t218 * (pkin(2) * t364 - t254 * t256 - t254) + t276 * t252) * t246 - t103 * t365;
t306 = t79 * t399;
t317 = pkin(2) * t365;
t401 = t136 * t94;
t82 = t256 * t401;
t414 = t266 * t401 + t82;
t417 = -(-t414 * t365 + (t79 * t317 + (t422 * t79 + (t265 + (1 + t256) * pkin(2)) * t94) * t246) * t136) * t399 - (t246 * t422 + t317) * t136 * t306;
t31 = t300 * t136 + t417;
t427 = -t221 * g(1) - t218 * g(2);
t28 = t31 + t427;
t16 = -t170 * t246 + t28 * t249;
t241 = legFrame(2,3);
t219 = sin(t241);
t222 = cos(t241);
t171 = g(1) * t219 - g(2) * t222;
t247 = sin(qJ(2,2));
t250 = cos(qJ(2,2));
t257 = (qJ(3,2) ^ 2);
t213 = -t257 - t422;
t236 = t250 ^ 2;
t338 = t247 * t250;
t314 = pkin(2) * t338;
t146 = 0.2e1 * qJ(3,2) * t314 + t213 * t236 - t257 + t422;
t139 = 0.1e1 / t146;
t374 = qJ(3,2) * t219;
t195 = pkin(2) * t374;
t151 = t213 * t222 - 0.2e1 * t195;
t373 = qJ(3,2) * t222;
t311 = pkin(2) * t373;
t154 = t213 * t219 + 0.2e1 * t311;
t107 = t151 * t236 + t154 * t338 - t222 * t266 + t195 - t222;
t110 = -t151 * t338 + t154 * t236 - t219 * t266 - t219 - t311;
t260 = koppelP(2,2);
t263 = koppelP(2,1);
t165 = t227 * t263 + t228 * t260;
t168 = -t227 * t260 + t228 * t263;
t125 = -t165 * t238 + t168 * t243 + t244;
t128 = -t165 * t243 - t168 * t238 + t245;
t299 = t107 * t128 + t110 * t125;
t140 = 0.1e1 / t146 ^ 2;
t369 = qJ(3,2) * t260;
t190 = pkin(2) * t263 + t369;
t368 = qJ(3,2) * t263;
t191 = pkin(2) * t260 - t368;
t113 = (t190 * t228 - t191 * t227) * t222 + (t190 * t227 + t191 * t228) * t219;
t370 = qJ(3,2) * t254;
t371 = qJ(3,2) * t253;
t104 = (t226 - t370) * t222 + t219 * (-t371 - t421) + t113 * t252;
t290 = t165 * t219 + t168 * t222;
t372 = qJ(3,2) * t250;
t326 = -0.2e1 * t372;
t95 = t104 * t247 + (-t219 * t254 + t253 * t222 + t290 * t252) * t326;
t394 = t140 * t95;
t230 = t257 * t260;
t318 = pkin(2) * t368;
t175 = t230 + t260 + t318;
t216 = pkin(2) * t369;
t335 = t257 * t263;
t176 = t216 - t263 - t335;
t275 = (t175 * t227 + t176 * t228) * t222 - t219 * (t175 * t228 - t176 * t227);
t80 = ((-pkin(2) * t370 - t253 * t257 - t253) * t222 - t219 * (pkin(2) * t371 - t254 * t257 - t254) + t275 * t252) * t247 - t104 * t372;
t305 = t80 * t394;
t319 = pkin(2) * t372;
t396 = t139 * t95;
t83 = t257 * t396;
t413 = t266 * t396 + t83;
t416 = -(-t413 * t372 + (t80 * t319 + (t422 * t80 + (t265 + (1 + t257) * pkin(2)) * t95) * t247) * t139) * t394 - (t247 * t422 + t319) * t139 * t305;
t32 = t299 * t139 + t416;
t426 = -t222 * g(1) - t219 * g(2);
t29 = t32 + t426;
t17 = -t171 * t247 + t29 * t250;
t242 = legFrame(1,3);
t220 = sin(t242);
t223 = cos(t242);
t172 = g(1) * t220 - g(2) * t223;
t248 = sin(qJ(2,1));
t251 = cos(qJ(2,1));
t258 = (qJ(3,1) ^ 2);
t214 = -t258 - t422;
t237 = t251 ^ 2;
t337 = t248 * t251;
t313 = pkin(2) * t337;
t147 = 0.2e1 * qJ(3,1) * t313 + t214 * t237 - t258 + t422;
t142 = 0.1e1 / t147;
t381 = qJ(3,1) * t220;
t196 = pkin(2) * t381;
t152 = t214 * t223 - 0.2e1 * t196;
t380 = qJ(3,1) * t223;
t312 = pkin(2) * t380;
t155 = t214 * t220 + 0.2e1 * t312;
t108 = t152 * t237 + t155 * t337 - t223 * t266 + t196 - t223;
t111 = -t152 * t337 + t155 * t237 - t220 * t266 - t220 - t312;
t261 = koppelP(1,2);
t264 = koppelP(1,1);
t166 = t227 * t264 + t228 * t261;
t169 = -t227 * t261 + t228 * t264;
t126 = -t166 * t238 + t169 * t243 + t244;
t129 = -t166 * t243 - t169 * t238 + t245;
t298 = t108 * t129 + t111 * t126;
t143 = 0.1e1 / t147 ^ 2;
t376 = qJ(3,1) * t261;
t192 = pkin(2) * t264 + t376;
t375 = qJ(3,1) * t264;
t193 = pkin(2) * t261 - t375;
t114 = (t192 * t228 - t193 * t227) * t223 + (t192 * t227 + t193 * t228) * t220;
t377 = qJ(3,1) * t254;
t378 = qJ(3,1) * t253;
t105 = (t226 - t377) * t223 + t220 * (-t378 - t421) + t114 * t252;
t289 = t166 * t220 + t169 * t223;
t379 = qJ(3,1) * t251;
t327 = -0.2e1 * t379;
t96 = t105 * t248 + (-t220 * t254 + t253 * t223 + t289 * t252) * t327;
t389 = t143 * t96;
t231 = t258 * t261;
t320 = pkin(2) * t375;
t177 = t231 + t261 + t320;
t217 = pkin(2) * t376;
t334 = t258 * t264;
t178 = t217 - t264 - t334;
t274 = (t177 * t227 + t178 * t228) * t223 - t220 * (t177 * t228 - t178 * t227);
t81 = ((-pkin(2) * t377 - t253 * t258 - t253) * t223 - t220 * (pkin(2) * t378 - t254 * t258 - t254) + t274 * t252) * t248 - t105 * t379;
t304 = t81 * t389;
t321 = pkin(2) * t379;
t391 = t142 * t96;
t84 = t258 * t391;
t412 = t266 * t391 + t84;
t415 = -(-t412 * t379 + (t81 * t321 + (t422 * t81 + (t265 + (1 + t258) * pkin(2)) * t96) * t248) * t142) * t389 - (t248 * t422 + t321) * t142 * t304;
t33 = t298 * t142 + t415;
t425 = -t223 * g(1) - t220 * g(2);
t30 = t33 + t425;
t18 = -t172 * t248 + t30 * t251;
t285 = t220 * t258 - t312;
t331 = t258 * t223 + t196;
t117 = t331 * t251 - t248 * (-t220 - t285);
t120 = t285 * t251 - t248 * (t223 + t331);
t295 = t117 * t129 + t120 * t126;
t430 = t295 * t142;
t284 = t219 * t257 - t311;
t332 = t257 * t222 + t195;
t116 = t332 * t250 - t247 * (-t219 - t284);
t119 = t284 * t250 - t247 * (t222 + t332);
t296 = t116 * t128 + t119 * t125;
t429 = t296 * t139;
t283 = t218 * t256 - t310;
t333 = t256 * t221 + t194;
t115 = t333 * t249 - t246 * (-t218 - t283);
t118 = t283 * t249 - t246 * (t221 + t333);
t297 = t115 * t127 + t118 * t124;
t428 = t297 * t136;
t21 = -t172 * t251 - t248 * t30;
t20 = -t171 * t250 - t247 * t29;
t19 = -t170 * t249 - t246 * t28;
t424 = pkin(2) * g(1);
t423 = pkin(2) * g(2);
t91 = t94 ^ 2;
t400 = t137 * t91;
t138 = t136 * t137;
t309 = t138 * t79 * t94;
t131 = t221 * t325 + t246 * (pkin(2) * t221 - t367);
t350 = t131 * t136;
t130 = 0.2e1 * t218 * t365 - t246 * (pkin(2) * t218 + t366);
t351 = t130 * t136;
t384 = t79 * t136;
t89 = pkin(2) * t401;
t68 = t89 - t384;
t49 = t127 * t351 + t124 * t350 - (-(pkin(2) * t68 - t82) * t339 + (-t384 + (0.2e1 * t89 - t384) * t235) * qJ(3,3)) * t399 - (-qJ(3,3) * t235 - qJ(3,3) + t315) * t309;
t43 = t246 * t49 + t249 * t400;
t330 = t422 * t262;
t158 = 0.2e1 * t215 - t330 - t336;
t303 = t422 * t259;
t159 = -t229 - t303 - 0.2e1 * t316;
t121 = t158 * t228 - t159 * t227;
t179 = t215 - t330;
t180 = -t303 - t316;
t294 = t158 * t227 + t159 * t228;
t76 = (t121 * t218 - t294 * t221) * t235 - (t121 * t221 + t294 * t218) * t339 + (t227 * t179 + t180 * t228) * t221 - (t179 * t228 - t227 * t180) * t218;
t420 = t43 * t76;
t92 = t95 ^ 2;
t395 = t140 * t92;
t141 = t139 * t140;
t308 = t141 * t80 * t95;
t133 = t222 * t326 + t247 * (pkin(2) * t222 - t374);
t348 = t133 * t139;
t132 = 0.2e1 * t219 * t372 - t247 * (pkin(2) * t219 + t373);
t349 = t132 * t139;
t383 = t80 * t139;
t90 = pkin(2) * t396;
t69 = t90 - t383;
t50 = t128 * t349 + t125 * t348 - (-(pkin(2) * t69 - t83) * t338 + (-t383 + (0.2e1 * t90 - t383) * t236) * qJ(3,2)) * t394 - (-qJ(3,2) * t236 - qJ(3,2) + t314) * t308;
t44 = t247 * t50 + t250 * t395;
t329 = t422 * t263;
t160 = 0.2e1 * t216 - t329 - t335;
t302 = t422 * t260;
t161 = -t230 - t302 - 0.2e1 * t318;
t122 = t160 * t228 - t161 * t227;
t181 = t216 - t329;
t182 = -t302 - t318;
t293 = t160 * t227 + t161 * t228;
t77 = (t122 * t219 - t293 * t222) * t236 - (t122 * t222 + t293 * t219) * t338 + (t227 * t181 + t182 * t228) * t222 - (t181 * t228 - t227 * t182) * t219;
t419 = t44 * t77;
t93 = t96 ^ 2;
t390 = t143 * t93;
t144 = t142 * t143;
t307 = t144 * t81 * t96;
t135 = t223 * t327 + t248 * (pkin(2) * t223 - t381);
t346 = t135 * t142;
t134 = 0.2e1 * t220 * t379 - t248 * (pkin(2) * t220 + t380);
t347 = t134 * t142;
t382 = t81 * t142;
t85 = pkin(2) * t391;
t67 = t85 - t382;
t51 = t129 * t347 + t126 * t346 - (-(pkin(2) * t67 - t84) * t337 + (-t382 + (0.2e1 * t85 - t382) * t237) * qJ(3,1)) * t389 - (-qJ(3,1) * t237 - qJ(3,1) + t313) * t307;
t45 = t248 * t51 + t251 * t390;
t328 = t422 * t264;
t162 = 0.2e1 * t217 - t328 - t334;
t301 = t422 * t261;
t163 = -t231 - t301 - 0.2e1 * t320;
t123 = t162 * t228 - t163 * t227;
t183 = t217 - t328;
t184 = -t301 - t320;
t292 = t162 * t227 + t163 * t228;
t78 = (t123 * t220 - t292 * t223) * t237 - (t123 * t223 + t292 * t220) * t337 + (t227 * t183 + t184 * t228) * t223 - (t183 * t228 - t227 * t184) * t220;
t418 = t45 * t78;
t411 = qJ(3,1) * t51;
t410 = qJ(3,2) * t50;
t409 = qJ(3,3) * t49;
t408 = t106 * t43;
t407 = t107 * t44;
t406 = t108 * t45;
t405 = t109 * t43;
t404 = t110 * t44;
t403 = t111 * t45;
t402 = t136 * t76;
t398 = t138 * t91;
t397 = t139 * t77;
t393 = t141 * t92;
t392 = t142 * t78;
t388 = t144 * t93;
t100 = t112 * t246 + t291 * t325;
t360 = t100 * t136;
t101 = t113 * t247 + t290 * t326;
t359 = t101 * t139;
t102 = t114 * t248 + t289 * t327;
t358 = t102 * t142;
t357 = t106 * t136;
t356 = t107 * t139;
t355 = t108 * t142;
t354 = t109 * t136;
t353 = t110 * t139;
t352 = t111 * t142;
t40 = -t246 * t400 + t249 * t49;
t41 = -t247 * t395 + t250 * t50;
t42 = -t248 * t390 + t251 * t51;
t324 = t42 * t392 + t41 * t397 + t40 * t402;
t323 = t42 * t355 + t41 * t356 + t40 * t357;
t322 = t42 * t352 + t41 * t353 + t40 * t354;
t70 = 0.2e1 * t306;
t71 = 0.2e1 * t305;
t72 = 0.2e1 * t304;
t288 = -(pkin(2) * qJ(3,3) + t339) * t309 + (t68 * t339 + ((-pkin(2) * t79 - t235 * t94 + t94) * t136 + t414) * qJ(3,3)) * t399;
t287 = -(pkin(2) * qJ(3,2) + t338) * t308 + (t69 * t338 + ((-pkin(2) * t80 - t236 * t95 + t95) * t139 + t413) * qJ(3,2)) * t394;
t286 = -(pkin(2) * qJ(3,1) + t337) * t307 + (t67 * t337 + ((-pkin(2) * t81 - t237 * t96 + t96) * t142 + t412) * qJ(3,1)) * t389;
t48 = t51 * pkin(2);
t282 = -qJ(3,1) * t390 - t286 - t48;
t47 = t50 * pkin(2);
t281 = -qJ(3,2) * t395 - t287 - t47;
t46 = t49 * pkin(2);
t280 = -qJ(3,3) * t400 - t288 - t46;
t279 = t288 - t428;
t278 = t287 - t429;
t277 = t286 - t430;
t225 = t245 - g(1);
t224 = t244 - g(2);
t205 = -g(1) * qJ(3,1) + t423;
t204 = -g(2) * qJ(3,1) - t424;
t203 = -g(1) * qJ(3,2) + t423;
t202 = -g(2) * qJ(3,2) - t424;
t201 = -g(1) * qJ(3,3) + t423;
t200 = -g(2) * qJ(3,3) - t424;
t157 = -t227 * t243 - t228 * t238;
t156 = -t227 * t238 + t228 * t243;
t149 = t224 * t227 + t225 * t228;
t148 = t224 * t228 - t225 * t227;
t99 = -t114 * t379 + t274 * t248;
t98 = -t113 * t372 + t275 * t247;
t97 = -t112 * t365 + t276 * t246;
t15 = t72 + 0.2e1 * t411 - t21;
t14 = t71 + 0.2e1 * t410 - t20;
t13 = t70 + 0.2e1 * t409 - t19;
t12 = t18 + t277 + 0.2e1 * t48;
t11 = t17 + t278 + 0.2e1 * t47;
t10 = t16 + t279 + 0.2e1 * t46;
t9 = t282 + t430 - t18;
t8 = t281 + t429 - t17;
t7 = t280 + t428 - t16;
t6 = -t282 * t251 + (-pkin(2) * t390 + t411 + t72) * t248 + (-t295 * t251 + t298) * t142 + t415 + t425;
t5 = -t281 * t250 + (-pkin(2) * t395 + t410 + t71) * t247 + (-t296 * t250 + t299) * t139 + t416 + t426;
t4 = -t280 * t249 + (-pkin(2) * t400 + t409 + t70) * t246 + (-t297 * t249 + t300) * t136 + t417 + t427;
t3 = (pkin(2) * t33 + t204 * t223 - t205 * t220) * t251 + (qJ(3,1) * t33 + t204 * t220 + t205 * t223) * t248 + t258 * t51 + qJ(3,1) * t72 + pkin(2) * (t48 + t277);
t2 = (pkin(2) * t32 + t202 * t222 - t203 * t219) * t250 + (qJ(3,2) * t32 + t202 * t219 + t203 * t222) * t247 + t257 * t50 + qJ(3,2) * t71 + pkin(2) * (t47 + t278);
t1 = (pkin(2) * t31 + t200 * t221 - t201 * t218) * t249 + (qJ(3,3) * t31 + t200 * t218 + t201 * t221) * t246 + t256 * t49 + qJ(3,3) * t70 + pkin(2) * (t46 + t279);
t22 = [t28 * t357 + t29 * t356 + t30 * t355, t51 * t347 + t349 * t50 + t351 * t49, t16 * t351 + t17 * t349 + t18 * t347 + t323 (t134 * t21 - t406) * t142 + (t132 * t20 - t407) * t139 + (t130 * t19 - t408) * t136 (-t117 * t51 + t12 * t134) * t142 + (t11 * t132 - t116 * t50) * t139 + (t10 * t130 - t115 * t49) * t136 + t323, -t115 * t398 - t116 * t393 - t117 * t388 + (t134 * t15 + t406) * t142 + (t132 * t14 + t407) * t139 + (t13 * t130 + t408) * t136 (t108 * t6 + t117 * t9 + t134 * t3) * t142 + (t107 * t5 + t116 * t8 + t132 * t2) * t139 + (t1 * t130 + t106 * t4 + t115 * t7) * t136, 0, t157, -t156, -t148 * t227 + t149 * t228; t28 * t354 + t29 * t353 + t30 * t352, t346 * t51 + t348 * t50 + t350 * t49, t16 * t350 + t17 * t348 + t18 * t346 + t322 (t135 * t21 - t403) * t142 + (t133 * t20 - t404) * t139 + (t131 * t19 - t405) * t136 (t12 * t135 - t120 * t51) * t142 + (t11 * t133 - t119 * t50) * t139 + (t10 * t131 - t118 * t49) * t136 + t322, -t118 * t398 - t119 * t393 - t120 * t388 + (t135 * t15 + t403) * t142 + (t133 * t14 + t404) * t139 + (t13 * t131 + t405) * t136 (t111 * t6 + t120 * t9 + t135 * t3) * t142 + (t110 * t5 + t119 * t8 + t133 * t2) * t139 + (t1 * t131 + t109 * t4 + t118 * t7) * t136, 0, t156, t157, t148 * t228 + t149 * t227; t28 * t402 + t29 * t397 + t30 * t392, t358 * t51 + t359 * t50 + t360 * t49, t16 * t360 + t17 * t359 + t18 * t358 + t324 (t102 * t21 - t418) * t142 + (t101 * t20 - t419) * t139 + (t100 * t19 - t420) * t136 (t102 * t12 - t51 * t99) * t142 + (t101 * t11 - t50 * t98) * t139 + (t10 * t100 - t49 * t97) * t136 + t324, -t97 * t398 - t98 * t393 - t99 * t388 + (t102 * t15 + t418) * t142 + (t101 * t14 + t419) * t139 + (t100 * t13 + t420) * t136 (t102 * t3 + t6 * t78 + t9 * t99) * t142 + (t101 * t2 + t5 * t77 + t8 * t98) * t139 + (t1 * t100 + t4 * t76 + t7 * t97) * t136, t243, t148, -t149, 0;];
tauX_reg  = t22;
