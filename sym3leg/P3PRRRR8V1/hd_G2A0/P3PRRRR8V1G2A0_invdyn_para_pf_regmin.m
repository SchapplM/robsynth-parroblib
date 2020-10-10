% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V1G2A0
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
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR8V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:03:02
% EndTime: 2020-08-06 17:03:11
% DurationCPUTime: 8.10s
% Computational Cost: add. (33936->402), mult. (88266->837), div. (6471->14), fcn. (86574->22), ass. (0->360)
t200 = cos(qJ(3,3));
t211 = 0.1e1 / pkin(2);
t201 = cos(qJ(2,3));
t195 = sin(qJ(2,3));
t326 = t195 * t200;
t152 = pkin(2) * t326 - t201 * pkin(5);
t185 = sin(pkin(3));
t187 = cos(pkin(3));
t194 = sin(qJ(3,3));
t391 = pkin(2) * t194;
t277 = t152 * t185 + t187 * t391;
t395 = 0.1e1 / t277;
t349 = t395 * t211;
t295 = t395 * t349;
t188 = legFrame(3,2);
t172 = sin(t188);
t175 = cos(t188);
t207 = xDP(2);
t208 = xDP(1);
t147 = -t172 * t207 + t175 * t208;
t186 = cos(pkin(6));
t206 = xDP(3);
t163 = t206 * t186;
t184 = sin(pkin(6));
t116 = t184 * t147 + t163;
t330 = t187 * t206;
t336 = t187 * t195;
t340 = t185 * t200;
t85 = ((-t147 * t336 - t201 * t206) * t184 + (t201 * t147 - t195 * t330) * t186) * t194 - t116 * t340;
t271 = t85 * t295;
t341 = t184 * t206;
t113 = t147 * t186 - t341;
t333 = t187 * t201;
t388 = pkin(2) * t200;
t79 = -(t113 * t195 + t116 * t333) * t388 + (t113 * t201 - t116 * t336) * pkin(5);
t250 = t79 * t271;
t178 = 0.1e1 / t200;
t213 = t200 ^ 2;
t179 = 0.1e1 / t213;
t191 = xDDP(3);
t192 = xDDP(2);
t193 = xDDP(1);
t291 = t178 * t349;
t243 = t185 * t79 * t291;
t328 = t194 * t195;
t221 = t187 * t328 + t340;
t327 = t194 * t201;
t106 = t221 * t184 - t186 * t327;
t351 = t395 * t178;
t293 = t175 * t351;
t258 = t106 * t293;
t294 = t172 * t351;
t259 = t106 * t294;
t109 = -t184 * t327 - t221 * t186;
t299 = t109 * t351;
t372 = t395 * t85;
t311 = t395 * t372;
t319 = pkin(5) * t372;
t329 = t187 * t211;
t350 = t395 * t194;
t123 = 0.1e1 / t277 ^ 2;
t381 = t123 * t79;
t382 = t395 * t79;
t280 = t194 * t319;
t67 = (t280 - t382) * t178;
t25 = t191 * t299 + t192 * t259 - t193 * t258 + (-(t187 * t67 + (pkin(2) * (t185 * t201 * t372 + t329 * t382) * t213 - t185 * (t79 * t350 - t319) * t326) * t178) * t311 + (-t201 * t243 + (t185 * t328 + (t178 - t200) * t187) * t372) * t381) * t179;
t357 = t25 * t194;
t394 = t179 - 0.2e1;
t13 = t200 * t357 - t394 * t250;
t400 = 0.2e1 * t13;
t202 = cos(qJ(3,2));
t203 = cos(qJ(2,2));
t197 = sin(qJ(2,2));
t323 = t197 * t202;
t153 = pkin(2) * t323 - t203 * pkin(5);
t196 = sin(qJ(3,2));
t390 = pkin(2) * t196;
t276 = t153 * t185 + t187 * t390;
t396 = 0.1e1 / t276;
t346 = t396 * t211;
t290 = t396 * t346;
t189 = legFrame(2,2);
t173 = sin(t189);
t176 = cos(t189);
t148 = -t173 * t207 + t176 * t208;
t117 = t184 * t148 + t163;
t335 = t187 * t197;
t339 = t185 * t202;
t86 = ((-t148 * t335 - t203 * t206) * t184 + (t203 * t148 - t197 * t330) * t186) * t196 - t117 * t339;
t269 = t86 * t290;
t114 = t148 * t186 - t341;
t332 = t187 * t203;
t387 = pkin(2) * t202;
t80 = -(t114 * t197 + t117 * t332) * t387 + (t114 * t203 - t117 * t335) * pkin(5);
t249 = t80 * t269;
t180 = 0.1e1 / t202;
t214 = t202 ^ 2;
t181 = 0.1e1 / t214;
t286 = t180 * t346;
t241 = t185 * t80 * t286;
t325 = t196 * t197;
t220 = t187 * t325 + t339;
t324 = t196 * t203;
t107 = t220 * t184 - t186 * t324;
t348 = t396 * t180;
t288 = t176 * t348;
t256 = t107 * t288;
t289 = t173 * t348;
t257 = t107 * t289;
t110 = -t184 * t324 - t220 * t186;
t298 = t110 * t348;
t368 = t396 * t86;
t309 = t396 * t368;
t318 = pkin(5) * t368;
t347 = t396 * t196;
t125 = 0.1e1 / t276 ^ 2;
t377 = t125 * t80;
t378 = t396 * t80;
t279 = t196 * t318;
t68 = (t279 - t378) * t180;
t26 = t191 * t298 + t192 * t257 - t193 * t256 + (-(t187 * t68 + (pkin(2) * (t185 * t203 * t368 + t329 * t378) * t214 - t185 * (t80 * t347 - t318) * t323) * t180) * t309 + (-t203 * t241 + (t185 * t325 + (t180 - t202) * t187) * t368) * t377) * t181;
t356 = t26 * t196;
t393 = t181 - 0.2e1;
t14 = t202 * t356 - t393 * t249;
t399 = 0.2e1 * t14;
t204 = cos(qJ(3,1));
t205 = cos(qJ(2,1));
t199 = sin(qJ(2,1));
t320 = t199 * t204;
t154 = pkin(2) * t320 - t205 * pkin(5);
t198 = sin(qJ(3,1));
t389 = pkin(2) * t198;
t275 = t154 * t185 + t187 * t389;
t397 = 0.1e1 / t275;
t343 = t397 * t211;
t285 = t397 * t343;
t190 = legFrame(1,2);
t174 = sin(t190);
t177 = cos(t190);
t149 = -t174 * t207 + t177 * t208;
t118 = t184 * t149 + t163;
t334 = t187 * t199;
t338 = t185 * t204;
t87 = ((-t149 * t334 - t205 * t206) * t184 + (t205 * t149 - t199 * t330) * t186) * t198 - t118 * t338;
t267 = t87 * t285;
t115 = t149 * t186 - t341;
t331 = t187 * t205;
t386 = pkin(2) * t204;
t81 = -(t115 * t199 + t118 * t331) * t386 + (t115 * t205 - t118 * t334) * pkin(5);
t248 = t81 * t267;
t322 = t198 * t199;
t219 = t187 * t322 + t338;
t321 = t198 * t205;
t108 = t219 * t184 - t186 * t321;
t182 = 0.1e1 / t204;
t215 = t204 ^ 2;
t183 = 0.1e1 / t215;
t281 = t182 * t343;
t239 = t185 * t81 * t281;
t345 = t397 * t182;
t283 = t177 * t345;
t255 = t108 * t283;
t342 = t174 * t182;
t284 = t397 * t342;
t111 = -t184 * t321 - t219 * t186;
t297 = t111 * t345;
t364 = t397 * t87;
t307 = t397 * t364;
t317 = pkin(5) * t364;
t344 = t397 * t198;
t127 = 0.1e1 / t275 ^ 2;
t373 = t127 * t81;
t374 = t397 * t81;
t278 = t198 * t317;
t69 = (t278 - t374) * t182;
t27 = t108 * t192 * t284 + t191 * t297 - t193 * t255 + (-(t187 * t69 + (pkin(2) * (t185 * t205 * t364 + t329 * t374) * t215 - t185 * (t81 * t344 - t317) * t320) * t182) * t307 + (-t205 * t239 + (t185 * t322 + (t182 - t204) * t187) * t364) * t373) * t183;
t355 = t27 * t198;
t392 = t183 - 0.2e1;
t15 = t204 * t355 - t392 * t248;
t398 = 0.2e1 * t15;
t385 = t184 * g(3);
t384 = t395 * t25;
t171 = t186 * g(3);
t253 = g(1) * t175 - g(2) * t172;
t143 = t253 * t184 + t171;
t337 = t186 * t187;
t150 = -t185 * g(1) + g(2) * t337;
t151 = g(1) * t337 + t185 * g(2);
t159 = t187 * t385;
t209 = pkin(5) ^ 2;
t210 = pkin(2) ^ 2;
t292 = t178 * t350;
t64 = -pkin(5) * t79 * t292 + (t178 * t209 + t200 * t210) * t372;
t155 = pkin(5) * t195 + t201 * t388;
t224 = -t152 * t187 + t185 * t391;
t94 = -t184 * t155 + t224 * t186;
t88 = t94 * t172 + t175 * t277;
t89 = t172 * t277 - t94 * t175;
t97 = t186 * t155 + t224 * t184;
t55 = (t191 * t97 + t192 * t88 + t193 * t89) * t395 + (t64 * t311 - t67 * t381) * t178;
t218 = t150 * t172 - t151 * t175 + t55 * t185 + t159;
t43 = t195 * t143 + t218 * t201;
t383 = t395 * t43;
t380 = t396 * t26;
t252 = g(1) * t176 - g(2) * t173;
t144 = t252 * t184 + t171;
t287 = t180 * t347;
t65 = -pkin(5) * t80 * t287 + (t180 * t209 + t202 * t210) * t368;
t156 = pkin(5) * t197 + t203 * t387;
t223 = -t153 * t187 + t185 * t390;
t95 = -t184 * t156 + t223 * t186;
t90 = t95 * t173 + t176 * t276;
t91 = t173 * t276 - t95 * t176;
t98 = t186 * t156 + t223 * t184;
t56 = (t191 * t98 + t192 * t90 + t193 * t91) * t396 + (t65 * t309 - t68 * t377) * t180;
t217 = t150 * t173 - t151 * t176 + t56 * t185 + t159;
t44 = t197 * t144 + t217 * t203;
t379 = t396 * t44;
t376 = t397 * t27;
t251 = g(1) * t177 - g(2) * t174;
t145 = t251 * t184 + t171;
t282 = t182 * t344;
t66 = -pkin(5) * t81 * t282 + (t182 * t209 + t204 * t210) * t364;
t157 = pkin(5) * t199 + t205 * t386;
t222 = -t154 * t187 + t185 * t389;
t96 = -t184 * t157 + t222 * t186;
t92 = t96 * t174 + t177 * t275;
t93 = t174 * t275 - t96 * t177;
t99 = t186 * t157 + t222 * t184;
t57 = (t191 * t99 + t192 * t92 + t193 * t93) * t397 + (t66 * t307 - t69 * t373) * t182;
t216 = t150 * t174 - t151 * t177 + t57 * t185 + t159;
t45 = t199 * t145 + t216 * t205;
t375 = t397 * t45;
t371 = t395 * t88;
t370 = t395 * t89;
t369 = t395 * t97;
t367 = t396 * t90;
t366 = t396 * t91;
t365 = t396 * t98;
t363 = t397 * t92;
t362 = t397 * t93;
t361 = t397 * t99;
t360 = t201 * t25;
t359 = t203 * t26;
t358 = t205 * t27;
t354 = t85 ^ 2 * t123;
t353 = t86 ^ 2 * t125;
t352 = t87 ^ 2 * t127;
t103 = (t184 * t333 + t186 * t195) * t388 + (t184 * t336 - t186 * t201) * pkin(5);
t316 = t103 * t384;
t104 = (t184 * t332 + t186 * t197) * t387 + (t184 * t335 - t186 * t203) * pkin(5);
t315 = t104 * t380;
t146 = t184 * t334 - t186 * t205;
t105 = (t184 * t331 + t186 * t199) * t386 + t146 * pkin(5);
t314 = t105 * t376;
t313 = t106 * t383;
t312 = t107 * t379;
t212 = 0.1e1 / pkin(2) ^ 2;
t310 = t123 * t212 * t79 ^ 2;
t308 = t125 * t212 * t80 ^ 2;
t306 = t127 * t212 * t81 ^ 2;
t305 = t179 * t354;
t304 = t181 * t353;
t303 = t183 * t352;
t100 = -(-t184 * t195 + t186 * t333) * t388 - pkin(5) * (t184 * t201 + t186 * t336);
t302 = t100 * t351;
t101 = -(-t184 * t197 + t186 * t332) * t387 - pkin(5) * (t184 * t203 + t186 * t335);
t301 = t101 * t348;
t102 = -(-t184 * t199 + t186 * t331) * t386 - pkin(5) * (t184 * t205 + t186 * t334);
t300 = t102 * t345;
t112 = t146 * t198 + t184 * t338;
t296 = t112 * t342;
t274 = t43 * t299;
t273 = t44 * t298;
t272 = t45 * t297;
t270 = t25 * t292;
t268 = t26 * t287;
t266 = t27 * t282;
t265 = t103 * t294;
t264 = t103 * t293;
t263 = t104 * t289;
t262 = t104 * t288;
t261 = t105 * t284;
t260 = t105 * t283;
t254 = t397 * t296;
t247 = t103 * t270;
t246 = t104 * t268;
t245 = t105 * t266;
t244 = t305 * t350;
t242 = t304 * t347;
t240 = t303 * t344;
t238 = t43 * t259;
t237 = t44 * t257;
t236 = t43 * t258;
t235 = t44 * t256;
t234 = t45 * t255;
t233 = (t310 + t354) * t179 * t195 - t360;
t232 = (t308 + t353) * t181 * t197 - t359;
t231 = (t306 + t352) * t183 * t199 - t358;
t230 = 0.2e1 * t250;
t229 = 0.2e1 * t249;
t228 = 0.2e1 * t248;
t227 = t103 * t244;
t226 = t104 * t242;
t225 = t105 * t240;
t142 = t251 * t186 - t385;
t141 = t252 * t186 - t385;
t140 = t253 * t186 - t385;
t121 = t205 * t145;
t120 = t203 * t144;
t119 = t201 * t143;
t75 = t392 * t352;
t74 = t393 * t353;
t73 = t394 * t354;
t54 = -t174 * g(1) - t177 * g(2) + t57;
t53 = -t173 * g(1) - t176 * g(2) + t56;
t52 = -t172 * g(1) - t175 * g(2) + t55;
t51 = -t185 * t142 - t54 * t187;
t50 = -t185 * t141 - t53 * t187;
t49 = -t185 * t140 - t52 * t187;
t48 = t121 + (t142 * t187 - t185 * t54) * t199;
t47 = t120 + (t141 * t187 - t185 * t53) * t197;
t46 = t119 + (t140 * t187 - t185 * t52) * t195;
t42 = -t216 * t199 + t121;
t41 = -t217 * t197 + t120;
t40 = -t218 * t195 + t119;
t39 = (-t187 * t66 * t267 - (-t198 * t154 * t239 + t187 * (-t182 * t278 + t204 * t374)) * t81 * t285) * t183 + (t102 * t191 + (t174 * t192 - t177 * t193) * t105) * t281;
t38 = (-t187 * t65 * t269 - (-t196 * t153 * t241 + t187 * (-t180 * t279 + t202 * t378)) * t80 * t290) * t181 + (t101 * t191 + (t173 * t192 - t176 * t193) * t104) * t286;
t37 = (-t187 * t64 * t271 - (-t194 * t152 * t243 + t187 * (-t178 * t280 + t200 * t382)) * t79 * t295) * t179 + (t100 * t191 + (t172 * t192 - t175 * t193) * t103) * t291;
t36 = t182 * t306 + t39 * t198;
t35 = t180 * t308 + t38 * t196;
t34 = t178 * t310 + t37 * t194;
t33 = -t183 * t198 * t306 + t39 * t204;
t32 = -t181 * t196 * t308 + t38 * t202;
t31 = -t179 * t194 * t310 + t37 * t200;
t30 = t183 * t205 * t228 + t199 * t39;
t29 = t181 * t203 * t229 + t197 * t38;
t28 = t179 * t201 * t230 + t195 * t37;
t24 = t197 * t26 + t203 * t304;
t23 = t197 * t304 - t359;
t22 = t199 * t27 + t205 * t303;
t21 = t195 * t25 + t201 * t305;
t20 = t199 * t303 - t358;
t19 = t195 * t305 - t360;
t18 = (t182 * t228 + t355) * t198;
t17 = (t180 * t229 + t356) * t196;
t16 = (t178 * t230 + t357) * t194;
t12 = t51 * t198 + t48 * t204;
t11 = t48 * t198 - t51 * t204;
t10 = t50 * t196 + t47 * t202;
t9 = t47 * t196 - t50 * t202;
t8 = t49 * t194 + t46 * t200;
t7 = t46 * t194 - t49 * t200;
t6 = (t231 * t198 - t204 * t30) * t185 - t187 * t36;
t5 = (t232 * t196 - t202 * t29) * t185 - t187 * t35;
t4 = (t233 * t194 - t200 * t28) * t185 - t187 * t34;
t3 = (-t198 * t30 - t231 * t204) * t185 + t187 * t33;
t2 = (-t196 * t29 - t232 * t202) * t185 + t187 * t32;
t1 = (-t194 * t28 - t233 * t200) * t185 + t187 * t31;
t58 = [t54 * t362 + t53 * t366 + t52 * t370, -t25 * t258 - t27 * t255 - t26 * t256, -t236 - t235 - t234 + (-t19 * t370 - t20 * t362 - t23 * t366) * t185, -t40 * t258 - t41 * t256 - t42 * t255 + (-t21 * t370 - t22 * t362 - t24 * t366) * t185, -t16 * t258 - t17 * t256 - t18 * t255 + (t175 * t227 + t176 * t226 + t177 * t225) * t211, (-t260 * t75 - t262 * t74 - t264 * t73) * t211 - 0.2e1 * t13 * t258 - 0.2e1 * t14 * t256 - 0.2e1 * t15 * t255, -t34 * t258 - t35 * t256 - t36 * t255 + (-t175 * t247 - t176 * t246 - t177 * t245) * t211, -t31 * t258 - t32 * t256 - t33 * t255 + (-t175 * t316 - t176 * t315 - t177 * t314) * t211, (-t260 * t39 - t262 * t38 - t264 * t37) * t211, -t175 * t313 - t176 * t312 - t177 * t108 * t375 + t1 * t370 + t2 * t366 + t3 * t362 + (-t11 * t260 - t262 * t9 - t264 * t7) * t211, t194 * t236 + t196 * t235 + t198 * t234 + t4 * t370 + t5 * t366 + t6 * t362 + (-t10 * t262 - t12 * t260 - t264 * t8) * t211, t193 - g(1); t54 * t363 + t53 * t367 + t52 * t371, t25 * t259 + t27 * t254 + t26 * t257, t238 + t237 + t45 * t254 + (-t19 * t371 - t20 * t363 - t23 * t367) * t185, t40 * t259 + t41 * t257 + t42 * t254 + (-t21 * t371 - t22 * t363 - t24 * t367) * t185, t16 * t259 + t17 * t257 + t18 * t254 + (-t172 * t227 - t173 * t226 - t174 * t225) * t211, (t261 * t75 + t263 * t74 + t265 * t73) * t211 + t259 * t400 + t257 * t399 + t254 * t398, t34 * t259 + t35 * t257 + t36 * t254 + (t172 * t247 + t173 * t246 + t174 * t245) * t211, t31 * t259 + t32 * t257 + t33 * t254 + (t172 * t316 + t173 * t315 + t174 * t314) * t211, (t261 * t39 + t263 * t38 + t265 * t37) * t211, t172 * t313 + t173 * t312 + t1 * t371 + t2 * t367 + (t112 * t174 * t45 + t92 * t3) * t397 + (t11 * t261 + t263 * t9 + t265 * t7) * t211, -t194 * t238 - t196 * t237 + t4 * t371 + t5 * t367 + (-t198 * t296 * t45 + t92 * t6) * t397 + (t10 * t263 + t12 * t261 + t265 * t8) * t211, t192 - g(2); t54 * t361 + t53 * t365 + t52 * t369, t25 * t299 + t26 * t298 + t27 * t297, t274 + t273 + t272 + (-t19 * t369 - t20 * t361 - t23 * t365) * t185, t40 * t299 + t41 * t298 + t42 * t297 + (-t21 * t369 - t22 * t361 - t24 * t365) * t185, t16 * t299 + t17 * t298 + t18 * t297 + (-t100 * t244 - t101 * t242 - t102 * t240) * t211, (t300 * t75 + t301 * t74 + t302 * t73) * t211 + t299 * t400 + t298 * t399 + t297 * t398, t34 * t299 + t35 * t298 + t36 * t297 + (t100 * t270 + t101 * t268 + t102 * t266) * t211, t31 * t299 + t32 * t298 + t33 * t297 + (t100 * t384 + t101 * t380 + t102 * t376) * t211, (t300 * t39 + t301 * t38 + t302 * t37) * t211, t1 * t369 + t109 * t383 + t110 * t379 + t111 * t375 + t2 * t365 + t3 * t361 + (t11 * t300 + t301 * t9 + t302 * t7) * t211, -t194 * t274 - t196 * t273 - t198 * t272 + t4 * t369 + t5 * t365 + t6 * t361 + (t10 * t301 + t12 * t300 + t302 * t8) * t211, t191 - g(3);];
tauX_reg  = t58;
