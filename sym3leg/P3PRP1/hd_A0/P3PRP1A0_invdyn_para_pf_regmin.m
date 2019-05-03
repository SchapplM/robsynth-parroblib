% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRP1A0
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
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRP1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:42:01
% EndTime: 2019-05-03 14:42:10
% DurationCPUTime: 9.50s
% Computational Cost: add. (78519->521), mult. (156556->778), div. (2844->6), fcn. (50847->14), ass. (0->369)
t259 = (pkin(2) ^ 2);
t429 = -t259 - 1;
t232 = legFrame(3,3);
t212 = sin(t232);
t215 = cos(t232);
t173 = t215 * g(1) + t212 * g(2);
t238 = sin(qJ(2,3));
t241 = cos(qJ(2,3));
t249 = (qJ(3,3) ^ 2);
t206 = -t249 - t429;
t227 = t241 ^ 2;
t340 = t241 * t238;
t308 = pkin(2) * t340;
t151 = 0.2e1 * qJ(3,3) * t308 + t206 * t227 - t249 + t429;
t142 = 0.1e1 / t151;
t372 = qJ(3,3) * t212;
t316 = pkin(2) * t372;
t156 = t206 * t215 + 0.2e1 * t316;
t371 = qJ(3,3) * t215;
t311 = pkin(2) * t371;
t278 = -t212 * t206 + 0.2e1 * t311;
t103 = t156 * t227 - t215 * t259 + t278 * t340 - t215 - t316;
t104 = -t156 * t340 + t259 * t212 + t278 * t227 + t212 - t311;
t248 = xP(3);
t222 = sin(t248);
t223 = cos(t248);
t252 = koppelP(3,2);
t255 = koppelP(3,1);
t167 = t222 * t255 + t223 * t252;
t170 = -t222 * t252 + t223 * t255;
t244 = xDP(3);
t230 = t244 ^ 2;
t235 = xDDP(3);
t236 = xDDP(2);
t127 = -t230 * t167 + t170 * t235 + t236;
t237 = xDDP(1);
t130 = -t167 * t235 - t230 * t170 + t237;
t302 = t103 * t127 + t104 * t130;
t258 = pkin(2) * t259;
t143 = 0.1e1 / t151 ^ 2;
t245 = xDP(2);
t246 = xDP(1);
t293 = t167 * t215 - t212 * t170;
t367 = qJ(3,3) * t252;
t191 = pkin(2) * t255 - t367;
t366 = qJ(3,3) * t255;
t192 = pkin(2) * t252 + t366;
t139 = t191 * t223 - t222 * t192;
t220 = t245 * pkin(2);
t221 = pkin(2) * t246;
t290 = t222 * t191 + t192 * t223;
t368 = qJ(3,3) * t246;
t369 = qJ(3,3) * t245;
t296 = t212 * (t139 * t244 + t220 + t368) - (t290 * t244 - t221 + t369) * t215;
t370 = qJ(3,3) * t241;
t327 = 0.2e1 * t370;
t94 = t296 * t238 + (-t245 * t212 - t246 * t215 + t293 * t244) * t327;
t408 = t143 * t94;
t319 = pkin(2) * t366;
t338 = t249 * t252;
t176 = -t252 + t319 - t338;
t209 = pkin(2) * t367;
t177 = t249 * t255 + t209 + t255;
t269 = (-t176 * t223 + t222 * t177) * t215 - (t222 * t176 + t177 * t223) * t212;
t79 = ((-pkin(2) * t369 - t249 * t246 - t246) * t215 - (-pkin(2) * t368 + t249 * t245 + t245) * t212 + t269 * t244) * t238 - t296 * t370;
t305 = t79 * t408;
t315 = pkin(2) * t370;
t388 = t94 * t142;
t82 = t249 * t388;
t419 = t259 * t388 + t82;
t422 = -(-t419 * t370 + (t79 * t315 + t238 * (t429 * t79 + (t258 + (1 + t249) * pkin(2)) * t94)) * t142) * t408 - (t238 * t429 + t315) * t142 * t305;
t31 = t302 * t142 + t422;
t434 = t212 * g(1) - t215 * g(2);
t28 = t31 + t434;
t16 = t173 * t238 + t28 * t241;
t233 = legFrame(2,3);
t213 = sin(t233);
t216 = cos(t233);
t174 = t216 * g(1) + t213 * g(2);
t239 = sin(qJ(2,2));
t242 = cos(qJ(2,2));
t250 = (qJ(3,2) ^ 2);
t207 = -t250 - t429;
t228 = t242 ^ 2;
t339 = t242 * t239;
t307 = pkin(2) * t339;
t152 = 0.2e1 * qJ(3,2) * t307 + t207 * t228 - t250 + t429;
t145 = 0.1e1 / t152;
t379 = qJ(3,2) * t213;
t317 = pkin(2) * t379;
t157 = t207 * t216 + 0.2e1 * t317;
t378 = qJ(3,2) * t216;
t310 = pkin(2) * t378;
t281 = -t213 * t207 + 0.2e1 * t310;
t105 = t157 * t228 - t216 * t259 + t281 * t339 - t216 - t317;
t106 = -t157 * t339 + t259 * t213 + t281 * t228 + t213 - t310;
t253 = koppelP(2,2);
t256 = koppelP(2,1);
t168 = t222 * t256 + t223 * t253;
t171 = -t222 * t253 + t223 * t256;
t128 = -t230 * t168 + t171 * t235 + t236;
t131 = -t168 * t235 - t230 * t171 + t237;
t301 = t105 * t128 + t106 * t131;
t146 = 0.1e1 / t152 ^ 2;
t292 = t168 * t216 - t213 * t171;
t374 = qJ(3,2) * t253;
t193 = pkin(2) * t256 - t374;
t373 = qJ(3,2) * t256;
t194 = pkin(2) * t253 + t373;
t140 = t193 * t223 - t222 * t194;
t289 = t222 * t193 + t194 * t223;
t375 = qJ(3,2) * t246;
t376 = qJ(3,2) * t245;
t295 = t213 * (t140 * t244 + t220 + t375) - (t289 * t244 - t221 + t376) * t216;
t377 = qJ(3,2) * t242;
t329 = 0.2e1 * t377;
t95 = t295 * t239 + (-t245 * t213 - t246 * t216 + t292 * t244) * t329;
t404 = t146 * t95;
t320 = pkin(2) * t373;
t337 = t250 * t253;
t178 = -t253 + t320 - t337;
t210 = pkin(2) * t374;
t179 = t250 * t256 + t210 + t256;
t268 = (-t178 * t223 + t222 * t179) * t216 - (t222 * t178 + t179 * t223) * t213;
t80 = ((-pkin(2) * t376 - t250 * t246 - t246) * t216 - (-pkin(2) * t375 + t250 * t245 + t245) * t213 + t268 * t244) * t239 - t295 * t377;
t304 = t80 * t404;
t321 = pkin(2) * t377;
t387 = t95 * t145;
t83 = t250 * t387;
t418 = t259 * t387 + t83;
t421 = -(-t418 * t377 + (t80 * t321 + t239 * (t429 * t80 + (t258 + (1 + t250) * pkin(2)) * t95)) * t145) * t404 - (t239 * t429 + t321) * t145 * t304;
t32 = t301 * t145 + t421;
t435 = t213 * g(1) - t216 * g(2);
t29 = t32 + t435;
t17 = t174 * t239 + t29 * t242;
t234 = legFrame(1,3);
t214 = sin(t234);
t217 = cos(t234);
t175 = t217 * g(1) + t214 * g(2);
t240 = sin(qJ(2,1));
t243 = cos(qJ(2,1));
t251 = (qJ(3,1) ^ 2);
t208 = -t251 - t429;
t229 = t243 ^ 2;
t341 = t240 * t243;
t306 = pkin(2) * t341;
t153 = 0.2e1 * qJ(3,1) * t306 + t208 * t229 - t251 + t429;
t148 = 0.1e1 / t153;
t386 = qJ(3,1) * t214;
t318 = pkin(2) * t386;
t158 = t208 * t217 + 0.2e1 * t318;
t385 = qJ(3,1) * t217;
t309 = pkin(2) * t385;
t284 = -t214 * t208 + 0.2e1 * t309;
t107 = t158 * t229 - t217 * t259 + t284 * t341 - t217 - t318;
t108 = -t158 * t341 + t259 * t214 + t284 * t229 + t214 - t309;
t254 = koppelP(1,2);
t257 = koppelP(1,1);
t169 = t222 * t257 + t223 * t254;
t172 = -t222 * t254 + t223 * t257;
t129 = -t230 * t169 + t172 * t235 + t236;
t132 = -t169 * t235 - t230 * t172 + t237;
t300 = t107 * t129 + t108 * t132;
t149 = 0.1e1 / t153 ^ 2;
t291 = t169 * t217 - t214 * t172;
t381 = qJ(3,1) * t254;
t195 = pkin(2) * t257 - t381;
t380 = qJ(3,1) * t257;
t196 = pkin(2) * t254 + t380;
t141 = t195 * t223 - t222 * t196;
t288 = t222 * t195 + t196 * t223;
t382 = qJ(3,1) * t246;
t383 = qJ(3,1) * t245;
t294 = t214 * (t141 * t244 + t220 + t382) - (t288 * t244 - t221 + t383) * t217;
t384 = qJ(3,1) * t243;
t331 = 0.2e1 * t384;
t96 = t294 * t240 + (-t245 * t214 - t246 * t217 + t291 * t244) * t331;
t399 = t149 * t96;
t322 = pkin(2) * t380;
t336 = t251 * t254;
t180 = -t254 + t322 - t336;
t211 = pkin(2) * t381;
t181 = t251 * t257 + t211 + t257;
t267 = (-t180 * t223 + t222 * t181) * t217 - (t222 * t180 + t181 * t223) * t214;
t81 = ((-pkin(2) * t383 - t251 * t246 - t246) * t217 - (-pkin(2) * t382 + t251 * t245 + t245) * t214 + t267 * t244) * t240 - t294 * t384;
t303 = t81 * t399;
t323 = pkin(2) * t384;
t401 = t148 * t96;
t84 = t251 * t401;
t417 = t259 * t401 + t84;
t420 = -(-t417 * t384 + (t81 * t323 + t240 * (t429 * t81 + (t258 + (1 + t251) * pkin(2)) * t96)) * t148) * t399 - (t240 * t429 + t323) * t148 * t303;
t33 = t300 * t148 + t420;
t436 = t214 * g(1) - t217 * g(2);
t30 = t33 + t436;
t18 = t175 * t240 + t30 * t243;
t282 = t217 * t251 - t318;
t283 = -t251 * t214 - t309;
t111 = t283 * t243 + t240 * (-t217 - t282);
t114 = t282 * t243 - t240 * (t214 - t283);
t297 = t111 * t132 + t114 * t129;
t439 = t297 * t148;
t279 = t216 * t250 - t317;
t280 = -t250 * t213 - t310;
t110 = t280 * t242 + t239 * (-t216 - t279);
t113 = t279 * t242 - t239 * (t213 - t280);
t298 = t110 * t131 + t113 * t128;
t438 = t298 * t145;
t276 = t215 * t249 - t316;
t277 = -t249 * t212 - t311;
t109 = t277 * t241 + t238 * (-t215 - t276);
t112 = t276 * t241 - t238 * (t212 - t277);
t299 = t109 * t130 + t112 * t127;
t437 = t299 * t142;
t21 = t175 * t243 - t30 * t240;
t20 = t174 * t242 - t29 * t239;
t19 = t173 * t241 - t28 * t238;
t433 = -t139 * t212 + t290 * t215;
t432 = -t140 * t213 + t289 * t216;
t431 = -t141 * t214 + t288 * t217;
t430 = pkin(2) * g(2);
t91 = t94 ^ 2;
t409 = t143 * t91;
t144 = t142 * t143;
t314 = t144 * t79 * t94;
t328 = -0.2e1 * t370;
t134 = t215 * t328 + t238 * (pkin(2) * t215 + t372);
t355 = t134 * t142;
t133 = t212 * t328 + t238 * (pkin(2) * t212 - t371);
t356 = t133 * t142;
t391 = t79 * t142;
t89 = pkin(2) * t388;
t68 = t89 - t391;
t49 = t130 * t355 + t127 * t356 - (-(pkin(2) * t68 - t82) * t340 + (-t391 + (0.2e1 * t89 - t391) * t227) * qJ(3,3)) * t408 - (-t227 * qJ(3,3) - qJ(3,3) + t308) * t314;
t43 = t49 * t238 + t241 * t409;
t161 = t206 * t255 - 0.2e1 * t209;
t335 = t429 * t252;
t162 = 0.2e1 * t319 - t335 - t338;
t121 = t161 * t223 - t222 * t162;
t122 = t222 * t161 + t162 * t223;
t182 = t259 * t255 - t209 + t255;
t183 = t319 - t335;
t76 = (t121 * t215 + t122 * t212) * t227 + (-t121 * t212 + t122 * t215) * t340 + (-t182 * t223 + t222 * t183) * t215 - (t222 * t182 + t183 * t223) * t212;
t425 = t43 * t76;
t92 = t95 ^ 2;
t405 = t146 * t92;
t147 = t145 * t146;
t313 = t147 * t80 * t95;
t330 = -0.2e1 * t377;
t136 = t216 * t330 + t239 * (pkin(2) * t216 + t379);
t353 = t136 * t145;
t135 = t213 * t330 + t239 * (pkin(2) * t213 - t378);
t354 = t135 * t145;
t390 = t80 * t145;
t90 = pkin(2) * t387;
t69 = t90 - t390;
t50 = t131 * t353 + t128 * t354 - (-(pkin(2) * t69 - t83) * t339 + (-t390 + (0.2e1 * t90 - t390) * t228) * qJ(3,2)) * t404 - (-t228 * qJ(3,2) - qJ(3,2) + t307) * t313;
t44 = t50 * t239 + t242 * t405;
t163 = t207 * t256 - 0.2e1 * t210;
t334 = t429 * t253;
t164 = 0.2e1 * t320 - t334 - t337;
t123 = t163 * t223 - t222 * t164;
t124 = t222 * t163 + t164 * t223;
t184 = t259 * t256 - t210 + t256;
t185 = t320 - t334;
t77 = (t123 * t216 + t124 * t213) * t228 + (-t123 * t213 + t124 * t216) * t339 + (-t184 * t223 + t222 * t185) * t216 - (t222 * t184 + t185 * t223) * t213;
t424 = t44 * t77;
t93 = t96 ^ 2;
t400 = t149 * t93;
t150 = t148 * t149;
t312 = t150 * t81 * t96;
t332 = -0.2e1 * t384;
t138 = t217 * t332 + t240 * (pkin(2) * t217 + t386);
t351 = t138 * t148;
t137 = t214 * t332 + t240 * (pkin(2) * t214 - t385);
t352 = t137 * t148;
t389 = t81 * t148;
t85 = pkin(2) * t401;
t67 = t85 - t389;
t51 = t132 * t351 + t129 * t352 - (-(pkin(2) * t67 - t84) * t341 + (-t389 + (0.2e1 * t85 - t389) * t229) * qJ(3,1)) * t399 - (-t229 * qJ(3,1) - qJ(3,1) + t306) * t312;
t45 = t51 * t240 + t243 * t400;
t165 = t208 * t257 - 0.2e1 * t211;
t333 = t429 * t254;
t166 = 0.2e1 * t322 - t333 - t336;
t125 = t165 * t223 - t222 * t166;
t126 = t222 * t165 + t166 * t223;
t186 = t259 * t257 - t211 + t257;
t187 = t322 - t333;
t78 = (t125 * t217 + t126 * t214) * t229 + (-t125 * t214 + t126 * t217) * t341 + (-t186 * t223 + t222 * t187) * t217 - (t222 * t186 + t187 * t223) * t214;
t423 = t45 * t78;
t416 = t103 * t43;
t415 = t104 * t43;
t414 = t105 * t44;
t413 = t106 * t44;
t412 = t107 * t45;
t411 = t108 * t45;
t410 = t142 * t76;
t407 = t144 * t91;
t406 = t145 * t77;
t403 = t147 * t92;
t402 = t148 * t78;
t398 = t150 * t93;
t394 = t49 * qJ(3,3);
t393 = t50 * qJ(3,2);
t392 = t51 * qJ(3,1);
t100 = -t433 * t238 + t293 * t327;
t365 = t100 * t142;
t101 = -t432 * t239 + t292 * t329;
t364 = t101 * t145;
t102 = -t431 * t240 + t291 * t331;
t363 = t102 * t148;
t362 = t103 * t142;
t361 = t104 * t142;
t360 = t105 * t145;
t359 = t106 * t145;
t358 = t107 * t148;
t357 = t108 * t148;
t40 = -t238 * t409 + t241 * t49;
t41 = -t239 * t405 + t242 * t50;
t42 = -t240 * t400 + t243 * t51;
t326 = t40 * t410 + t42 * t402 + t41 * t406;
t325 = t42 * t358 + t41 * t360 + t40 * t362;
t324 = t42 * t357 + t41 * t359 + t40 * t361;
t70 = 0.2e1 * t305;
t71 = 0.2e1 * t304;
t72 = 0.2e1 * t303;
t287 = -(qJ(3,3) * pkin(2) + t340) * t314 + (t68 * t340 + ((-pkin(2) * t79 - t227 * t94 + t94) * t142 + t419) * qJ(3,3)) * t408;
t286 = -(pkin(2) * qJ(3,2) + t339) * t313 + (t69 * t339 + ((-pkin(2) * t80 - t228 * t95 + t95) * t145 + t418) * qJ(3,2)) * t404;
t285 = -(pkin(2) * qJ(3,1) + t341) * t312 + (t67 * t341 + ((-pkin(2) * t81 - t229 * t96 + t96) * t148 + t417) * qJ(3,1)) * t399;
t48 = t51 * pkin(2);
t275 = -qJ(3,1) * t400 - t285 - t48;
t47 = t50 * pkin(2);
t274 = -qJ(3,2) * t405 - t286 - t47;
t46 = t49 * pkin(2);
t273 = -qJ(3,3) * t409 - t287 - t46;
t272 = t287 - t437;
t271 = t286 - t438;
t270 = t285 - t439;
t247 = pkin(2) * g(1);
t219 = t237 - g(1);
t218 = t236 - g(2);
t202 = g(1) * qJ(3,1) + t430;
t201 = -g(2) * qJ(3,1) + t247;
t200 = g(1) * qJ(3,2) + t430;
t199 = -g(2) * qJ(3,2) + t247;
t198 = g(1) * qJ(3,3) + t430;
t197 = -g(2) * qJ(3,3) + t247;
t160 = -t222 * t235 - t223 * t230;
t159 = -t222 * t230 + t223 * t235;
t155 = t222 * t218 + t223 * t219;
t154 = t223 * t218 - t222 * t219;
t99 = t267 * t240 + t431 * t384;
t98 = t268 * t239 + t432 * t377;
t97 = t269 * t238 + t433 * t370;
t15 = 0.2e1 * t392 + t72 - t21;
t14 = 0.2e1 * t393 + t71 - t20;
t13 = 0.2e1 * t394 + t70 - t19;
t12 = t18 + t270 + 0.2e1 * t48;
t11 = t17 + t271 + 0.2e1 * t47;
t10 = t16 + t272 + 0.2e1 * t46;
t9 = t275 + t439 - t18;
t8 = t274 + t438 - t17;
t7 = t273 + t437 - t16;
t6 = -t275 * t243 + (-pkin(2) * t400 + t392 + t72) * t240 + (-t297 * t243 + t300) * t148 + t420 + t436;
t5 = -t274 * t242 + (-pkin(2) * t405 + t393 + t71) * t239 + (-t298 * t242 + t301) * t145 + t421 + t435;
t4 = -t273 * t241 + (-pkin(2) * t409 + t394 + t70) * t238 + (-t299 * t241 + t302) * t142 + t422 + t434;
t3 = (pkin(2) * t33 + t201 * t214 - t202 * t217) * t243 + (t33 * qJ(3,1) + t201 * t217 + t202 * t214) * t240 + t251 * t51 + qJ(3,1) * t72 + pkin(2) * (t48 + t270);
t2 = (t32 * pkin(2) + t199 * t213 - t200 * t216) * t242 + (qJ(3,2) * t32 + t199 * t216 + t200 * t213) * t239 + t250 * t50 + qJ(3,2) * t71 + pkin(2) * (t47 + t271);
t1 = (t31 * pkin(2) + t197 * t212 - t198 * t215) * t241 + (qJ(3,3) * t31 + t197 * t215 + t198 * t212) * t238 + t249 * t49 + qJ(3,3) * t70 + pkin(2) * (t46 + t272);
t22 = [t28 * t361 + t29 * t359 + t30 * t357, t51 * t351 + t50 * t353 + t49 * t355, t16 * t355 + t17 * t353 + t18 * t351 + t324, (t138 * t21 - t411) * t148 + (t136 * t20 - t413) * t145 + (t134 * t19 - t415) * t142, (-t111 * t51 + t12 * t138) * t148 + (t11 * t136 - t110 * t50) * t145 + (t10 * t134 - t109 * t49) * t142 + t324, -t109 * t407 - t110 * t403 - t111 * t398 + (t138 * t15 + t411) * t148 + (t136 * t14 + t413) * t145 + (t13 * t134 + t415) * t142, (t108 * t6 + t111 * t9 + t138 * t3) * t148 + (t106 * t5 + t110 * t8 + t136 * t2) * t145 + (t1 * t134 + t104 * t4 + t109 * t7) * t142, 0, t160, -t159, -t222 * t154 + t223 * t155; t28 * t362 + t29 * t360 + t30 * t358, t51 * t352 + t50 * t354 + t49 * t356, t16 * t356 + t17 * t354 + t18 * t352 + t325, (t137 * t21 - t412) * t148 + (t135 * t20 - t414) * t145 + (t133 * t19 - t416) * t142, (-t114 * t51 + t12 * t137) * t148 + (t11 * t135 - t113 * t50) * t145 + (t10 * t133 - t112 * t49) * t142 + t325, -t112 * t407 - t113 * t403 - t114 * t398 + (t137 * t15 + t412) * t148 + (t135 * t14 + t414) * t145 + (t13 * t133 + t416) * t142, (t107 * t6 + t114 * t9 + t137 * t3) * t148 + (t105 * t5 + t113 * t8 + t135 * t2) * t145 + (t1 * t133 + t103 * t4 + t112 * t7) * t142, 0, t159, t160, t223 * t154 + t222 * t155; t28 * t410 + t29 * t406 + t30 * t402, t51 * t363 + t50 * t364 + t49 * t365, t16 * t365 + t17 * t364 + t18 * t363 + t326, (t102 * t21 - t423) * t148 + (t101 * t20 - t424) * t145 + (t100 * t19 - t425) * t142, (t102 * t12 - t51 * t99) * t148 + (t101 * t11 - t50 * t98) * t145 + (t10 * t100 - t49 * t97) * t142 + t326, -t97 * t407 - t98 * t403 - t99 * t398 + (t102 * t15 + t423) * t148 + (t101 * t14 + t424) * t145 + (t100 * t13 + t425) * t142, (t102 * t3 + t6 * t78 + t9 * t99) * t148 + (t101 * t2 + t5 * t77 + t8 * t98) * t145 + (t1 * t100 + t4 * t76 + t7 * t97) * t142, t235, t154, -t155, 0;];
tauX_reg  = t22;
