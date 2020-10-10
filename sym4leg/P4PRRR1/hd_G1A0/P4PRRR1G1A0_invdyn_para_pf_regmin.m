% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [4x11]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P4PRRR1G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:14:25
% EndTime: 2020-03-02 20:14:37
% DurationCPUTime: 12.38s
% Computational Cost: add. (47780->357), mult. (41343->628), div. (5856->11), fcn. (35904->42), ass. (0->310)
t240 = pkin(7) + qJ(2,1);
t223 = qJ(3,1) + t240;
t208 = sin(t223);
t211 = cos(t223);
t245 = legFrame(1,3);
t227 = sin(t245);
t231 = cos(t245);
t392 = -t208 * t227 + t231 * t211;
t239 = pkin(7) + qJ(2,2);
t222 = qJ(3,2) + t239;
t207 = sin(t222);
t210 = cos(t222);
t244 = legFrame(2,3);
t226 = sin(t244);
t230 = cos(t244);
t391 = -t207 * t226 + t230 * t210;
t238 = pkin(7) + qJ(2,3);
t221 = qJ(3,3) + t238;
t206 = sin(t221);
t209 = cos(t221);
t243 = legFrame(3,3);
t225 = sin(t243);
t229 = cos(t243);
t390 = -t206 * t225 + t229 * t209;
t237 = pkin(7) + qJ(2,4);
t214 = qJ(3,4) + t237;
t200 = sin(t214);
t201 = cos(t214);
t242 = legFrame(4,3);
t224 = sin(t242);
t228 = cos(t242);
t389 = -t200 * t224 + t228 * t201;
t271 = 0.1e1 / pkin(2);
t212 = sin(t237);
t213 = cos(t237);
t118 = t200 * t213 - t212 * t201;
t386 = 0.1e1 / t118;
t329 = t386 * t271;
t388 = -2 * pkin(2);
t387 = 2 * pkin(2);
t215 = sin(t238);
t218 = cos(t238);
t128 = t206 * t218 - t215 * t209;
t385 = 0.1e1 / t128;
t216 = sin(t239);
t219 = cos(t239);
t129 = t207 * t219 - t216 * t210;
t384 = 0.1e1 / t129;
t217 = sin(t240);
t220 = cos(t240);
t130 = t208 * t220 - t217 * t211;
t383 = 0.1e1 / t130;
t270 = 0.1e1 / pkin(3);
t260 = xP(4);
t235 = sin(t260);
t236 = cos(t260);
t261 = koppelP(4,2);
t265 = koppelP(4,1);
t155 = t235 * t265 + t236 * t261;
t257 = xDP(4);
t259 = xDP(1);
t139 = t155 * t257 - t259;
t276 = -pkin(2) * (t212 * t224 - t228 * t213) + t389 * pkin(3);
t336 = t276 * t139;
t159 = -t235 * t261 + t236 * t265;
t258 = xDP(2);
t143 = t159 * t257 + t258;
t132 = t228 * t200 + t224 * t201;
t76 = pkin(2) * (t228 * t212 + t224 * t213) + t132 * pkin(3);
t348 = t143 * t76;
t280 = (t336 - t348) * t270;
t55 = (-t139 * t228 + t224 * t143) * t201 + t200 * (t224 * t139 + t143 * t228);
t39 = (t55 + t280) * t329;
t46 = t280 * t329;
t303 = t271 * t39 * t46;
t289 = pkin(3) * t303;
t272 = 1 / pkin(2) ^ 2;
t340 = t55 * t386;
t302 = t272 * t340;
t284 = t200 * t212 + t201 * t213;
t33 = pkin(3) * t39 + t284 * t340;
t241 = t257 ^ 2;
t248 = xDDP(4);
t249 = xDDP(2);
t96 = -t241 * t155 + t159 * t248 + t249;
t68 = t132 * t96 * t329;
t250 = xDDP(1);
t100 = -t155 * t248 - t241 * t159 + t250;
t69 = t389 * t100 * t329;
t368 = t68 + t69;
t25 = -(-t33 * t302 - t289) * t386 + t368;
t382 = pkin(2) * t25;
t262 = koppelP(3,2);
t266 = koppelP(3,1);
t156 = t235 * t266 + t236 * t262;
t140 = t156 * t257 - t259;
t160 = -t235 * t262 + t236 * t266;
t144 = t160 * t257 + t258;
t59 = (-t140 * t229 + t225 * t144) * t209 + t206 * (t225 * t140 + t144 * t229);
t339 = t59 * t385;
t300 = t272 * t339;
t283 = t206 * t215 + t209 * t218;
t275 = -pkin(2) * (t215 * t225 - t229 * t218) + t390 * pkin(3);
t334 = t275 * t140;
t136 = t229 * t206 + t225 * t209;
t79 = pkin(2) * (t229 * t215 + t225 * t218) + t136 * pkin(3);
t347 = t144 * t79;
t279 = (t334 - t347) * t270;
t326 = t385 * t271;
t43 = (t59 + t279) * t326;
t34 = pkin(3) * t43 + t283 * t339;
t101 = -t156 * t248 - t241 * t160 + t250;
t97 = -t241 * t156 + t160 * t248 + t249;
t367 = (t101 * t390 + t136 * t97) * t326;
t47 = t279 * t326;
t374 = t43 * t47;
t378 = pkin(3) * t271;
t26 = -(-t34 * t300 - t374 * t378) * t385 + t367;
t381 = pkin(2) * t26;
t263 = koppelP(2,2);
t267 = koppelP(2,1);
t157 = t235 * t267 + t236 * t263;
t141 = t157 * t257 - t259;
t161 = -t235 * t263 + t236 * t267;
t145 = t161 * t257 + t258;
t60 = (-t141 * t230 + t226 * t145) * t210 + t207 * (t226 * t141 + t145 * t230);
t338 = t60 * t384;
t298 = t272 * t338;
t282 = t207 * t216 + t210 * t219;
t274 = -pkin(2) * (t216 * t226 - t230 * t219) + t391 * pkin(3);
t332 = t274 * t141;
t137 = t230 * t207 + t226 * t210;
t80 = pkin(2) * (t230 * t216 + t226 * t219) + t137 * pkin(3);
t346 = t145 * t80;
t278 = (t332 - t346) * t270;
t325 = t384 * t271;
t44 = (t60 + t278) * t325;
t35 = pkin(3) * t44 + t282 * t338;
t102 = -t157 * t248 - t241 * t161 + t250;
t98 = -t241 * t157 + t161 * t248 + t249;
t366 = (t102 * t391 + t137 * t98) * t325;
t48 = t278 * t325;
t373 = t44 * t48;
t27 = -(-t35 * t298 - t373 * t378) * t384 + t366;
t380 = pkin(2) * t27;
t264 = koppelP(1,2);
t268 = koppelP(1,1);
t158 = t235 * t268 + t236 * t264;
t142 = t158 * t257 - t259;
t162 = -t235 * t264 + t236 * t268;
t146 = t162 * t257 + t258;
t61 = (-t142 * t231 + t227 * t146) * t211 + t208 * (t227 * t142 + t146 * t231);
t337 = t61 * t383;
t296 = t272 * t337;
t281 = t208 * t217 + t211 * t220;
t273 = -pkin(2) * (t217 * t227 - t231 * t220) + t392 * pkin(3);
t330 = t273 * t142;
t138 = t231 * t208 + t227 * t211;
t81 = pkin(2) * (t231 * t217 + t227 * t220) + t138 * pkin(3);
t345 = t146 * t81;
t277 = (t330 - t345) * t270;
t324 = t383 * t271;
t45 = (t61 + t277) * t324;
t36 = t45 * pkin(3) + t281 * t337;
t103 = -t158 * t248 - t241 * t162 + t250;
t99 = -t241 * t158 + t162 * t248 + t249;
t365 = (t103 * t392 + t138 * t99) * t324;
t49 = t277 * t324;
t372 = t45 * t49;
t28 = -(-t36 * t296 - t372 * t378) * t383 + t365;
t379 = pkin(2) * t28;
t40 = (t59 + (t334 / 0.2e1 - t347 / 0.2e1) * t270) * t326;
t377 = t40 * t47;
t41 = (t60 + (t332 / 0.2e1 - t346 / 0.2e1) * t270) * t325;
t376 = t41 * t48;
t42 = (t61 + (t330 / 0.2e1 - t345 / 0.2e1) * t270) * t324;
t375 = t42 * t49;
t371 = t79 * t97;
t370 = t80 * t98;
t369 = t81 * t99;
t88 = t155 * t228 - t224 * t159;
t92 = t224 * t155 + t159 * t228;
t288 = -t92 * t200 + t88 * t201;
t364 = t386 * (pkin(2) * (-t92 * t212 + t88 * t213) + t288 * pkin(3));
t363 = t386 * t288;
t362 = t386 * t76;
t361 = t386 * t276;
t89 = t156 * t229 - t225 * t160;
t93 = t225 * t156 + t160 * t229;
t287 = -t93 * t206 + t89 * t209;
t360 = t385 * (pkin(2) * (-t93 * t215 + t89 * t218) + t287 * pkin(3));
t359 = t385 * t287;
t90 = t157 * t230 - t226 * t161;
t94 = t226 * t157 + t161 * t230;
t286 = -t94 * t207 + t90 * t210;
t358 = t384 * (pkin(2) * (-t94 * t216 + t90 * t219) + t286 * pkin(3));
t357 = t384 * t286;
t91 = t158 * t231 - t227 * t162;
t95 = t227 * t158 + t162 * t231;
t285 = -t95 * t208 + t91 * t211;
t356 = t383 * (pkin(2) * (-t95 * t217 + t91 * t220) + t285 * pkin(3));
t355 = t383 * t285;
t354 = t385 * t79;
t353 = t385 * t275;
t352 = t384 * t80;
t351 = t384 * t274;
t350 = t383 * t81;
t349 = t383 * t273;
t269 = pkin(3) ^ 2;
t305 = 0.2e1 * pkin(3);
t38 = (t55 + (t336 / 0.2e1 - t348 / 0.2e1) * t270) * t329;
t344 = t270 * (-t39 * t269 + (-t284 * t38 * t305 - t340) * pkin(2));
t343 = (-t43 * t269 + (-t283 * t40 * t305 - t339) * pkin(2)) * t270;
t342 = (-t44 * t269 + (-t282 * t41 * t305 - t338) * pkin(2)) * t270;
t341 = (-t45 * t269 + (-t281 * t42 * t305 - t337) * pkin(2)) * t270;
t335 = t275 * t101;
t333 = t274 * t102;
t331 = t273 * t103;
t328 = t386 * t389;
t327 = t386 * t132;
t323 = t385 * t390;
t322 = t385 * t136;
t321 = t384 * t391;
t320 = t384 * t137;
t319 = t383 * t392;
t318 = t383 * t138;
t196 = t242 + t214;
t188 = sin(t196);
t189 = cos(t196);
t309 = g(1) * t189 + g(2) * t188;
t197 = t243 + t221;
t190 = sin(t197);
t193 = cos(t197);
t308 = g(1) * t193 + g(2) * t190;
t198 = t244 + t222;
t191 = sin(t198);
t194 = cos(t198);
t307 = g(1) * t194 + g(2) * t191;
t199 = t245 + t223;
t192 = sin(t199);
t195 = cos(t199);
t306 = g(1) * t195 + g(2) * t192;
t304 = -(t284 * pkin(2) + pkin(3)) * t386 * t303 + (-t276 * t100 - t76 * t96) * t270 * t329;
t301 = 0.1e1 / t118 ^ 2 * t271 * t55 ^ 2;
t299 = 0.1e1 / t128 ^ 2 * t271 * t59 ^ 2;
t297 = 0.1e1 / t129 ^ 2 * t271 * t60 ^ 2;
t295 = 0.1e1 / t130 ^ 2 * t271 * t61 ^ 2;
t293 = g(1) * t188 - g(2) * t189;
t292 = g(1) * t190 - g(2) * t193;
t291 = g(1) * t191 - g(2) * t194;
t290 = g(1) * t192 - g(2) * t195;
t256 = cos(qJ(3,1));
t255 = cos(qJ(3,2));
t254 = cos(qJ(3,3));
t253 = sin(qJ(3,1));
t252 = sin(qJ(3,2));
t251 = sin(qJ(3,3));
t247 = cos(qJ(3,4));
t246 = sin(qJ(3,4));
t234 = t250 - g(1);
t233 = t249 - g(2);
t232 = (xDDP(3) - g(3));
t170 = t231 * g(1) + t227 * g(2);
t169 = t230 * g(1) + t226 * g(2);
t168 = t229 * g(1) + t225 * g(2);
t167 = t228 * g(1) + t224 * g(2);
t166 = t227 * g(1) - t231 * g(2);
t165 = t226 * g(1) - t230 * g(2);
t164 = t225 * g(1) - t229 * g(2);
t163 = t224 * g(1) - t228 * g(2);
t154 = -t235 * t248 - t236 * t241;
t153 = -t235 * t241 + t236 * t248;
t152 = t235 * t233 + t236 * t234;
t151 = t236 * t233 - t235 * t234;
t114 = t281 * pkin(2) + pkin(3);
t113 = t282 * pkin(2) + pkin(3);
t112 = t283 * pkin(2) + pkin(3);
t111 = -t166 * t217 + t170 * t220;
t110 = -t165 * t216 + t169 * t219;
t109 = -t164 * t215 + t168 * t218;
t108 = t166 * t220 + t170 * t217;
t107 = t165 * t219 + t169 * t216;
t106 = t164 * t218 + t168 * t215;
t105 = -t163 * t212 + t167 * t213;
t104 = t163 * t213 + t167 * t212;
t24 = -t253 * t379 + t256 * t295 + t306;
t23 = -t251 * t381 + t254 * t299 + t308;
t22 = t251 * t299 + t254 * t381 + t292;
t21 = t252 * t297 + t255 * t380 + t291;
t20 = -t252 * t380 + t255 * t297 + t307;
t19 = t253 * t295 + t256 * t379 + t290;
t18 = t246 * t301 + t247 * t382 + t293;
t17 = -t246 * t382 + t247 * t301 + t309;
t16 = -((-t36 - t341) * t296 + ((-pkin(3) + t114) * t372 + (t331 + t369) * t270) * t271) * t383 + t365;
t15 = -((-t35 - t342) * t298 + ((-pkin(3) + t113) * t373 + (t333 + t370) * t270) * t271) * t384 + t366;
t14 = -((-t34 - t343) * t300 + ((-pkin(3) + t112) * t374 + (t335 + t371) * t270) * t271) * t385 + t367;
t13 = -((-t36 - t341 / 0.2e1) * t296 + ((-pkin(3) + t114 / 0.2e1) * t372 + (t331 / 0.2e1 + t369 / 0.2e1) * t270) * t271) * t383 + t365;
t12 = -((-t35 - t342 / 0.2e1) * t298 + ((-pkin(3) + t113 / 0.2e1) * t373 + (t333 / 0.2e1 + t370 / 0.2e1) * t270) * t271) * t384 + t366;
t11 = -((-t34 - t343 / 0.2e1) * t300 + ((-pkin(3) + t112 / 0.2e1) * t374 + (t335 / 0.2e1 + t371 / 0.2e1) * t270) * t271) * t385 + t367;
t10 = -(-t289 + (-t33 - t344) * t302) * t386 + t304 + t368;
t9 = 0.2e1 * t68 + 0.2e1 * t69 - (-0.2e1 * t289 + (-0.2e1 * t33 - t344) * t302) * t386 + t304;
t8 = (t253 * t13 + t256 * t375) * t388 + t306;
t7 = (t252 * t12 + t255 * t376) * t388 + t307;
t6 = (t251 * t11 + t254 * t377) * t388 + t308;
t5 = (t13 * t256 - t253 * t375) * t387 + t290;
t4 = (t12 * t255 - t252 * t376) * t387 + t291;
t3 = (t11 * t254 - t251 * t377) * t387 + t292;
t2 = -pkin(2) * ((0.2e1 * t55 * t329 + t46) * t46 * t247 + t246 * t9) + t309;
t1 = pkin(2) * (-0.2e1 * t46 * t246 * t38 + t9 * t247) + t293;
t29 = [0, (t25 * t328 + t26 * t323 + t27 * t321 + t28 * t319) * t271, (t104 * t328 + t106 * t323 + t107 * t321 + t108 * t319) * t271, (t105 * t328 + t109 * t323 + t110 * t321 + t111 * t319) * t271, (t10 * t328 + t14 * t323 + t15 * t321 + t16 * t319 + (-t10 * t361 - t14 * t353 - t15 * t351 - t16 * t349) * t270) * t271, (t1 * t328 + t3 * t323 + t4 * t321 + t5 * t319 + (-t18 * t361 - t19 * t349 - t21 * t351 - t22 * t353) * t270) * t271, (t2 * t328 + t6 * t323 + t7 * t321 + t8 * t319 + (-t17 * t361 - t20 * t351 - t23 * t353 - t24 * t349) * t270) * t271, 0, t154, -t153, -t235 * t151 + t236 * t152; 0, (t25 * t327 + t26 * t322 + t27 * t320 + t28 * t318) * t271, (t104 * t327 + t106 * t322 + t107 * t320 + t108 * t318) * t271, (t105 * t327 + t109 * t322 + t110 * t320 + t111 * t318) * t271, (t10 * t327 + t14 * t322 + t15 * t320 + t16 * t318 + (-t10 * t362 - t14 * t354 - t15 * t352 - t16 * t350) * t270) * t271, (t1 * t327 + t3 * t322 + t4 * t320 + t5 * t318 + (-t18 * t362 - t19 * t350 - t21 * t352 - t22 * t354) * t270) * t271, (t2 * t327 + t6 * t322 + t7 * t320 + t8 * t318 + (-t17 * t362 - t20 * t352 - t23 * t354 - t24 * t350) * t270) * t271, 0, t153, t154, t236 * t151 + t235 * t152; 4 * t232, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232; 0, (-t25 * t363 - t26 * t359 - t27 * t357 - t28 * t355) * t271, (-t104 * t363 - t106 * t359 - t107 * t357 - t108 * t355) * t271, (-t105 * t363 - t109 * t359 - t110 * t357 - t111 * t355) * t271, (-t10 * t363 - t14 * t359 - t15 * t357 - t16 * t355 + (t10 * t364 + t14 * t360 + t15 * t358 + t16 * t356) * t270) * t271, (-t1 * t363 - t3 * t359 - t4 * t357 - t5 * t355 + (t18 * t364 + t19 * t356 + t21 * t358 + t22 * t360) * t270) * t271, (-t2 * t363 - t6 * t359 - t7 * t357 - t8 * t355 + (t17 * t364 + t20 * t358 + t23 * t360 + t24 * t356) * t270) * t271, t248, t151, -t152, 0;];
tauX_reg  = t29;
