% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRP1G1P1A0
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
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRP1G1P1A0_invdyn_para_pf_reg(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1G1P1A0_invdyn_para_pf_reg: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:35:07
% EndTime: 2018-12-20 17:35:17
% DurationCPUTime: 9.94s
% Computational Cost: add. (78519->521), mult. (156556->778), div. (2844->6), fcn. (50847->14), ass. (0->368)
t259 = (pkin(2) ^ 2);
t428 = -t259 - 1;
t232 = legFrame(3,3);
t212 = sin(t232);
t215 = cos(t232);
t173 = g(1) * t215 + g(2) * t212;
t238 = sin(qJ(2,3));
t241 = cos(qJ(2,3));
t249 = (qJ(3,3) ^ 2);
t206 = -t249 - t428;
t227 = t241 ^ 2;
t340 = t238 * t241;
t317 = pkin(2) * t340;
t151 = 0.2e1 * qJ(3,3) * t317 + t206 * t227 - t249 + t428;
t142 = 0.1e1 / t151;
t371 = qJ(3,3) * t212;
t312 = pkin(2) * t371;
t156 = t206 * t215 + 0.2e1 * t312;
t370 = qJ(3,3) * t215;
t308 = pkin(2) * t370;
t278 = -t206 * t212 + 0.2e1 * t308;
t103 = t156 * t227 - t215 * t259 + t278 * t340 - t215 - t312;
t104 = -t156 * t340 + t259 * t212 + t227 * t278 + t212 - t308;
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
t127 = -t167 * t230 + t170 * t235 + t236;
t237 = xDDP(1);
t130 = -t167 * t235 - t170 * t230 + t237;
t302 = t103 * t127 + t104 * t130;
t258 = pkin(2) * t259;
t143 = 0.1e1 / t151 ^ 2;
t245 = xDP(2);
t246 = xDP(1);
t293 = t167 * t215 - t170 * t212;
t366 = qJ(3,3) * t252;
t191 = pkin(2) * t255 - t366;
t365 = qJ(3,3) * t255;
t192 = pkin(2) * t252 + t365;
t139 = t191 * t223 - t192 * t222;
t220 = pkin(2) * t245;
t221 = pkin(2) * t246;
t290 = t191 * t222 + t192 * t223;
t367 = qJ(3,3) * t246;
t368 = qJ(3,3) * t245;
t296 = (t139 * t244 + t220 + t367) * t212 - (t244 * t290 - t221 + t368) * t215;
t369 = qJ(3,3) * t241;
t326 = 0.2e1 * t369;
t94 = t296 * t238 + (-t212 * t245 - t246 * t215 + t293 * t244) * t326;
t403 = t143 * t94;
t318 = pkin(2) * t365;
t337 = t249 * t252;
t176 = -t252 + t318 - t337;
t209 = pkin(2) * t366;
t177 = t249 * t255 + t209 + t255;
t269 = (-t176 * t223 + t177 * t222) * t215 - t212 * (t176 * t222 + t177 * t223);
t79 = ((-pkin(2) * t368 - t246 * t249 - t246) * t215 - t212 * (-pkin(2) * t367 + t245 * t249 + t245) + t269 * t244) * t238 - t296 * t369;
t305 = t79 * t403;
t319 = pkin(2) * t369;
t405 = t142 * t94;
t82 = t249 * t405;
t418 = t259 * t405 + t82;
t421 = -(-t418 * t369 + (t79 * t319 + (t428 * t79 + (t258 + (1 + t249) * pkin(2)) * t94) * t238) * t142) * t403 - (t238 * t428 + t319) * t142 * t305;
t31 = t302 * t142 + t421;
t433 = t212 * g(1) - t215 * g(2);
t28 = t31 + t433;
t16 = t173 * t238 + t28 * t241;
t233 = legFrame(2,3);
t213 = sin(t233);
t216 = cos(t233);
t174 = g(1) * t216 + g(2) * t213;
t239 = sin(qJ(2,2));
t242 = cos(qJ(2,2));
t250 = (qJ(3,2) ^ 2);
t207 = -t250 - t428;
t228 = t242 ^ 2;
t376 = qJ(3,2) * t242;
t313 = pkin(2) * t376;
t152 = t207 * t228 + 0.2e1 * t239 * t313 - t250 + t428;
t145 = 0.1e1 / t152;
t378 = qJ(3,2) * t213;
t314 = pkin(2) * t378;
t157 = t207 * t216 + 0.2e1 * t314;
t377 = qJ(3,2) * t216;
t307 = pkin(2) * t377;
t281 = -t207 * t213 + 0.2e1 * t307;
t339 = t239 * t242;
t105 = t157 * t228 - t216 * t259 + t281 * t339 - t216 - t314;
t106 = -t157 * t339 + t259 * t213 + t228 * t281 + t213 - t307;
t253 = koppelP(2,2);
t256 = koppelP(2,1);
t168 = t222 * t256 + t223 * t253;
t171 = -t222 * t253 + t223 * t256;
t128 = -t168 * t230 + t171 * t235 + t236;
t131 = -t168 * t235 - t171 * t230 + t237;
t301 = t105 * t128 + t106 * t131;
t146 = 0.1e1 / t152 ^ 2;
t292 = t168 * t216 - t171 * t213;
t373 = qJ(3,2) * t253;
t193 = pkin(2) * t256 - t373;
t372 = qJ(3,2) * t256;
t194 = pkin(2) * t253 + t372;
t140 = t193 * t223 - t194 * t222;
t289 = t193 * t222 + t194 * t223;
t374 = qJ(3,2) * t246;
t375 = qJ(3,2) * t245;
t295 = (t140 * t244 + t220 + t374) * t213 - (t244 * t289 - t221 + t375) * t216;
t328 = 0.2e1 * t376;
t95 = t295 * t239 + (-t213 * t245 - t246 * t216 + t292 * t244) * t328;
t398 = t146 * t95;
t320 = pkin(2) * t372;
t336 = t250 * t253;
t178 = -t253 + t320 - t336;
t210 = pkin(2) * t373;
t179 = t250 * t256 + t210 + t256;
t268 = (-t178 * t223 + t179 * t222) * t216 - t213 * (t178 * t222 + t179 * t223);
t80 = ((-pkin(2) * t375 - t246 * t250 - t246) * t216 - t213 * (-pkin(2) * t374 + t245 * t250 + t245) + t268 * t244) * t239 - t295 * t376;
t304 = t80 * t398;
t400 = t145 * t95;
t83 = t250 * t400;
t417 = t259 * t400 + t83;
t420 = -(-t417 * t376 + (t80 * t313 + (t428 * t80 + (t258 + (1 + t250) * pkin(2)) * t95) * t239) * t145) * t398 - (t239 * t428 + t313) * t145 * t304;
t32 = t301 * t145 + t420;
t434 = t213 * g(1) - t216 * g(2);
t29 = t32 + t434;
t17 = t174 * t239 + t29 * t242;
t234 = legFrame(1,3);
t214 = sin(t234);
t217 = cos(t234);
t175 = g(1) * t217 + g(2) * t214;
t240 = sin(qJ(2,1));
t243 = cos(qJ(2,1));
t251 = (qJ(3,1) ^ 2);
t208 = -t251 - t428;
t229 = t243 ^ 2;
t338 = t240 * t243;
t316 = pkin(2) * t338;
t153 = 0.2e1 * qJ(3,1) * t316 + t208 * t229 - t251 + t428;
t148 = 0.1e1 / t153;
t385 = qJ(3,1) * t214;
t315 = pkin(2) * t385;
t158 = t208 * t217 + 0.2e1 * t315;
t384 = qJ(3,1) * t217;
t306 = pkin(2) * t384;
t284 = -t208 * t214 + 0.2e1 * t306;
t107 = t158 * t229 - t217 * t259 + t284 * t338 - t217 - t315;
t108 = -t158 * t338 + t259 * t214 + t229 * t284 + t214 - t306;
t254 = koppelP(1,2);
t257 = koppelP(1,1);
t169 = t222 * t257 + t223 * t254;
t172 = -t222 * t254 + t223 * t257;
t129 = -t169 * t230 + t172 * t235 + t236;
t132 = -t169 * t235 - t172 * t230 + t237;
t300 = t107 * t129 + t108 * t132;
t149 = 0.1e1 / t153 ^ 2;
t291 = t169 * t217 - t172 * t214;
t380 = qJ(3,1) * t254;
t195 = pkin(2) * t257 - t380;
t379 = qJ(3,1) * t257;
t196 = pkin(2) * t254 + t379;
t141 = t195 * t223 - t196 * t222;
t288 = t195 * t222 + t196 * t223;
t381 = qJ(3,1) * t246;
t382 = qJ(3,1) * t245;
t294 = (t141 * t244 + t220 + t381) * t214 - (t244 * t288 - t221 + t382) * t217;
t383 = qJ(3,1) * t243;
t330 = 0.2e1 * t383;
t96 = t294 * t240 + (-t214 * t245 - t246 * t217 + t291 * t244) * t330;
t393 = t149 * t96;
t321 = pkin(2) * t379;
t335 = t251 * t254;
t180 = -t254 + t321 - t335;
t211 = pkin(2) * t380;
t181 = t251 * t257 + t211 + t257;
t267 = (-t180 * t223 + t181 * t222) * t217 - t214 * (t180 * t222 + t181 * t223);
t81 = ((-pkin(2) * t382 - t246 * t251 - t246) * t217 - t214 * (-pkin(2) * t381 + t245 * t251 + t245) + t267 * t244) * t240 - t294 * t383;
t303 = t81 * t393;
t322 = pkin(2) * t383;
t395 = t148 * t96;
t84 = t251 * t395;
t416 = t259 * t395 + t84;
t419 = -(-t416 * t383 + (t81 * t322 + (t428 * t81 + (t258 + (1 + t251) * pkin(2)) * t96) * t240) * t148) * t393 - (t240 * t428 + t322) * t148 * t303;
t33 = t300 * t148 + t419;
t435 = t214 * g(1) - t217 * g(2);
t30 = t33 + t435;
t18 = t175 * t240 + t30 * t243;
t282 = t217 * t251 - t315;
t283 = -t214 * t251 - t306;
t111 = t283 * t243 + t240 * (-t217 - t282);
t114 = t282 * t243 - t240 * (t214 - t283);
t297 = t111 * t132 + t114 * t129;
t438 = t297 * t148;
t279 = t216 * t250 - t314;
t280 = -t213 * t250 - t307;
t110 = t280 * t242 + t239 * (-t216 - t279);
t113 = t279 * t242 - t239 * (t213 - t280);
t298 = t110 * t131 + t113 * t128;
t437 = t298 * t145;
t276 = t215 * t249 - t312;
t277 = -t212 * t249 - t308;
t109 = t277 * t241 + t238 * (-t215 - t276);
t112 = t276 * t241 - t238 * (t212 - t277);
t299 = t109 * t130 + t112 * t127;
t436 = t299 * t142;
t21 = t175 * t243 - t240 * t30;
t20 = t174 * t242 - t239 * t29;
t19 = t173 * t241 - t238 * t28;
t432 = -t139 * t212 + t215 * t290;
t431 = -t140 * t213 + t216 * t289;
t430 = -t141 * t214 + t217 * t288;
t429 = pkin(2) * g(2);
t91 = t94 ^ 2;
t404 = t143 * t91;
t144 = t142 * t143;
t311 = t144 * t79 * t94;
t327 = -0.2e1 * t369;
t134 = t215 * t327 + t238 * (pkin(2) * t215 + t371);
t354 = t134 * t142;
t133 = t212 * t327 + t238 * (pkin(2) * t212 - t370);
t355 = t133 * t142;
t388 = t79 * t142;
t89 = pkin(2) * t405;
t68 = t89 - t388;
t49 = t130 * t354 + t127 * t355 - (-(pkin(2) * t68 - t82) * t340 + (-t388 + (0.2e1 * t89 - t388) * t227) * qJ(3,3)) * t403 - (-qJ(3,3) * t227 - qJ(3,3) + t317) * t311;
t43 = t49 * t238 + t241 * t404;
t161 = t206 * t255 - 0.2e1 * t209;
t334 = t428 * t252;
t162 = 0.2e1 * t318 - t334 - t337;
t121 = t161 * t223 - t162 * t222;
t122 = t161 * t222 + t162 * t223;
t182 = t255 * t259 - t209 + t255;
t183 = t318 - t334;
t76 = (t121 * t215 + t122 * t212) * t227 + (-t121 * t212 + t122 * t215) * t340 + (-t182 * t223 + t183 * t222) * t215 - (t182 * t222 + t183 * t223) * t212;
t424 = t43 * t76;
t92 = t95 ^ 2;
t399 = t146 * t92;
t147 = t145 * t146;
t310 = t147 * t80 * t95;
t329 = -0.2e1 * t376;
t136 = t216 * t329 + t239 * (pkin(2) * t216 + t378);
t352 = t136 * t145;
t135 = t213 * t329 + t239 * (pkin(2) * t213 - t377);
t353 = t135 * t145;
t387 = t80 * t145;
t90 = pkin(2) * t400;
t69 = t90 - t387;
t50 = t131 * t352 + t128 * t353 - (-(pkin(2) * t69 - t83) * t339 + (-t387 + (0.2e1 * t90 - t387) * t228) * qJ(3,2)) * t398 - (pkin(2) * t339 - qJ(3,2) * t228 - qJ(3,2)) * t310;
t44 = t50 * t239 + t242 * t399;
t163 = t207 * t256 - 0.2e1 * t210;
t333 = t428 * t253;
t164 = 0.2e1 * t320 - t333 - t336;
t123 = t163 * t223 - t164 * t222;
t124 = t163 * t222 + t164 * t223;
t184 = t256 * t259 - t210 + t256;
t185 = t320 - t333;
t77 = (t123 * t216 + t124 * t213) * t228 + (-t123 * t213 + t124 * t216) * t339 + (-t184 * t223 + t185 * t222) * t216 - (t184 * t222 + t185 * t223) * t213;
t423 = t44 * t77;
t93 = t96 ^ 2;
t394 = t149 * t93;
t150 = t148 * t149;
t309 = t150 * t81 * t96;
t331 = -0.2e1 * t383;
t138 = t217 * t331 + t240 * (pkin(2) * t217 + t385);
t350 = t138 * t148;
t137 = t214 * t331 + t240 * (pkin(2) * t214 - t384);
t351 = t137 * t148;
t386 = t81 * t148;
t85 = pkin(2) * t395;
t67 = t85 - t386;
t51 = t132 * t350 + t129 * t351 - (-(pkin(2) * t67 - t84) * t338 + (-t386 + (0.2e1 * t85 - t386) * t229) * qJ(3,1)) * t393 - (-qJ(3,1) * t229 - qJ(3,1) + t316) * t309;
t45 = t51 * t240 + t243 * t394;
t165 = t208 * t257 - 0.2e1 * t211;
t332 = t428 * t254;
t166 = 0.2e1 * t321 - t332 - t335;
t125 = t165 * t223 - t166 * t222;
t126 = t165 * t222 + t166 * t223;
t186 = t257 * t259 - t211 + t257;
t187 = t321 - t332;
t78 = (t125 * t217 + t126 * t214) * t229 + (-t125 * t214 + t126 * t217) * t338 + (-t186 * t223 + t187 * t222) * t217 - (t186 * t222 + t187 * t223) * t214;
t422 = t45 * t78;
t415 = qJ(3,1) * t51;
t414 = qJ(3,2) * t50;
t413 = qJ(3,3) * t49;
t412 = t103 * t43;
t411 = t104 * t43;
t410 = t105 * t44;
t409 = t106 * t44;
t408 = t107 * t45;
t407 = t108 * t45;
t406 = t142 * t76;
t402 = t144 * t91;
t401 = t145 * t77;
t397 = t147 * t92;
t396 = t148 * t78;
t392 = t150 * t93;
t100 = -t238 * t432 + t293 * t326;
t364 = t100 * t142;
t101 = -t239 * t431 + t292 * t328;
t363 = t101 * t145;
t102 = -t240 * t430 + t291 * t330;
t362 = t102 * t148;
t361 = t103 * t142;
t360 = t104 * t142;
t359 = t105 * t145;
t358 = t106 * t145;
t357 = t107 * t148;
t356 = t108 * t148;
t40 = -t238 * t404 + t241 * t49;
t41 = -t239 * t399 + t242 * t50;
t42 = -t240 * t394 + t243 * t51;
t325 = t42 * t396 + t40 * t406 + t41 * t401;
t324 = t42 * t357 + t41 * t359 + t40 * t361;
t323 = t42 * t356 + t41 * t358 + t40 * t360;
t70 = 0.2e1 * t305;
t71 = 0.2e1 * t304;
t72 = 0.2e1 * t303;
t287 = -(pkin(2) * qJ(3,3) + t340) * t311 + (t68 * t340 + ((-pkin(2) * t79 - t227 * t94 + t94) * t142 + t418) * qJ(3,3)) * t403;
t286 = -(qJ(3,2) * pkin(2) + t339) * t310 + (t69 * t339 + ((-pkin(2) * t80 - t228 * t95 + t95) * t145 + t417) * qJ(3,2)) * t398;
t285 = -(pkin(2) * qJ(3,1) + t338) * t309 + (t67 * t338 + ((-pkin(2) * t81 - t229 * t96 + t96) * t148 + t416) * qJ(3,1)) * t393;
t48 = t51 * pkin(2);
t275 = -qJ(3,1) * t394 - t285 - t48;
t47 = t50 * pkin(2);
t274 = -qJ(3,2) * t399 - t286 - t47;
t46 = t49 * pkin(2);
t273 = -qJ(3,3) * t404 - t287 - t46;
t272 = t287 - t436;
t271 = t286 - t437;
t270 = t285 - t438;
t247 = pkin(2) * g(1);
t219 = t237 - g(1);
t218 = t236 - g(2);
t202 = g(1) * qJ(3,1) + t429;
t201 = -g(2) * qJ(3,1) + t247;
t200 = g(1) * qJ(3,2) + t429;
t199 = -g(2) * qJ(3,2) + t247;
t198 = g(1) * qJ(3,3) + t429;
t197 = -g(2) * qJ(3,3) + t247;
t160 = -t222 * t235 - t223 * t230;
t159 = -t222 * t230 + t223 * t235;
t155 = t218 * t222 + t219 * t223;
t154 = t218 * t223 - t219 * t222;
t99 = t267 * t240 + t383 * t430;
t98 = t268 * t239 + t376 * t431;
t97 = t269 * t238 + t369 * t432;
t15 = t72 + 0.2e1 * t415 - t21;
t14 = t71 + 0.2e1 * t414 - t20;
t13 = t70 + 0.2e1 * t413 - t19;
t12 = t18 + t270 + 0.2e1 * t48;
t11 = t17 + t271 + 0.2e1 * t47;
t10 = t16 + t272 + 0.2e1 * t46;
t9 = t275 + t438 - t18;
t8 = t274 + t437 - t17;
t7 = t273 + t436 - t16;
t6 = -t275 * t243 + (-pkin(2) * t394 + t415 + t72) * t240 + (-t297 * t243 + t300) * t148 + t419 + t435;
t5 = -t274 * t242 + (-pkin(2) * t399 + t414 + t71) * t239 + (-t298 * t242 + t301) * t145 + t420 + t434;
t4 = -t273 * t241 + (-pkin(2) * t404 + t413 + t70) * t238 + (-t299 * t241 + t302) * t142 + t421 + t433;
t3 = (pkin(2) * t33 + t201 * t214 - t202 * t217) * t243 + (qJ(3,1) * t33 + t201 * t217 + t202 * t214) * t240 + t251 * t51 + qJ(3,1) * t72 + pkin(2) * (t48 + t270);
t2 = (pkin(2) * t32 + t199 * t213 - t200 * t216) * t242 + (qJ(3,2) * t32 + t199 * t216 + t200 * t213) * t239 + t250 * t50 + qJ(3,2) * t71 + pkin(2) * (t47 + t271);
t1 = (pkin(2) * t31 + t197 * t212 - t198 * t215) * t241 + (qJ(3,3) * t31 + t197 * t215 + t198 * t212) * t238 + t249 * t49 + qJ(3,3) * t70 + pkin(2) * (t46 + t272);
t22 = [t28 * t360 + t29 * t358 + t30 * t356, t51 * t350 + t50 * t352 + t49 * t354, t16 * t354 + t17 * t352 + t18 * t350 + t323 (t138 * t21 - t407) * t148 + (t136 * t20 - t409) * t145 + (t134 * t19 - t411) * t142 (-t111 * t51 + t12 * t138) * t148 + (t11 * t136 - t110 * t50) * t145 + (t10 * t134 - t109 * t49) * t142 + t323, -t109 * t402 - t110 * t397 - t111 * t392 + (t138 * t15 + t407) * t148 + (t136 * t14 + t409) * t145 + (t13 * t134 + t411) * t142 (t108 * t6 + t111 * t9 + t138 * t3) * t148 + (t106 * t5 + t110 * t8 + t136 * t2) * t145 + (t1 * t134 + t104 * t4 + t109 * t7) * t142, 0, t160, -t159, -t154 * t222 + t155 * t223; t28 * t361 + t29 * t359 + t30 * t357, t51 * t351 + t50 * t353 + t49 * t355, t16 * t355 + t17 * t353 + t18 * t351 + t324 (t137 * t21 - t408) * t148 + (t135 * t20 - t410) * t145 + (t133 * t19 - t412) * t142 (-t114 * t51 + t12 * t137) * t148 + (t11 * t135 - t113 * t50) * t145 + (t10 * t133 - t112 * t49) * t142 + t324, -t112 * t402 - t113 * t397 - t114 * t392 + (t137 * t15 + t408) * t148 + (t135 * t14 + t410) * t145 + (t13 * t133 + t412) * t142 (t107 * t6 + t114 * t9 + t137 * t3) * t148 + (t105 * t5 + t113 * t8 + t135 * t2) * t145 + (t1 * t133 + t103 * t4 + t112 * t7) * t142, 0, t159, t160, t154 * t223 + t155 * t222; t28 * t406 + t29 * t401 + t30 * t396, t51 * t362 + t50 * t363 + t49 * t364, t16 * t364 + t17 * t363 + t18 * t362 + t325 (t102 * t21 - t422) * t148 + (t101 * t20 - t423) * t145 + (t100 * t19 - t424) * t142 (t102 * t12 - t51 * t99) * t148 + (t101 * t11 - t50 * t98) * t145 + (t10 * t100 - t49 * t97) * t142 + t325, -t97 * t402 - t98 * t397 - t99 * t392 + (t102 * t15 + t422) * t148 + (t101 * t14 + t423) * t145 + (t100 * t13 + t424) * t142 (t102 * t3 + t6 * t78 + t9 * t99) * t148 + (t101 * t2 + t5 * t77 + t8 * t98) * t145 + (t1 * t100 + t4 * t76 + t7 * t97) * t142, t235, t154, -t155, 0;];
tauX_reg  = t22;
