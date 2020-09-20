% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V1G3A0
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR8V1G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:20
% EndTime: 2020-08-06 17:16:29
% DurationCPUTime: 8.90s
% Computational Cost: add. (33936->401), mult. (88305->827), div. (6471->14), fcn. (86599->22), ass. (0->356)
t198 = cos(qJ(3,3));
t209 = 0.1e1 / pkin(2);
t199 = cos(qJ(2,3));
t193 = sin(qJ(2,3));
t323 = t193 * t198;
t151 = pkin(2) * t323 - t199 * pkin(5);
t183 = sin(pkin(3));
t185 = cos(pkin(3));
t192 = sin(qJ(3,3));
t387 = pkin(2) * t192;
t275 = t151 * t183 + t185 * t387;
t393 = 0.1e1 / t275;
t345 = t393 * t209;
t293 = t393 * t345;
t186 = legFrame(3,2);
t170 = sin(t186);
t173 = cos(t186);
t205 = xDP(2);
t206 = xDP(1);
t146 = -t170 * t205 + t173 * t206;
t184 = cos(pkin(6));
t182 = sin(pkin(6));
t204 = xDP(3);
t337 = t182 * t204;
t113 = t146 * t184 - t337;
t327 = t185 * t204;
t333 = t185 * t193;
t336 = t183 * t198;
t85 = ((-t146 * t333 - t199 * t204) * t184 - (t199 * t146 - t193 * t327) * t182) * t192 - t113 * t336;
t269 = t85 * t293;
t161 = t204 * t184;
t116 = t182 * t146 + t161;
t330 = t185 * t199;
t384 = pkin(2) * t198;
t79 = -(t113 * t330 - t193 * t116) * t384 - (t113 * t333 + t199 * t116) * pkin(5);
t249 = t79 * t269;
t176 = 0.1e1 / t198;
t211 = t198 ^ 2;
t177 = 0.1e1 / t211;
t189 = xDDP(3);
t190 = xDDP(2);
t191 = xDDP(1);
t289 = t176 * t345;
t243 = t183 * t79 * t289;
t325 = t192 * t193;
t219 = t185 * t325 + t336;
t324 = t192 * t199;
t110 = t182 * t324 + t219 * t184;
t347 = t393 * t176;
t291 = t173 * t347;
t257 = t110 * t291;
t292 = t170 * t347;
t258 = t110 * t292;
t107 = t219 * t182 - t184 * t324;
t296 = t107 * t347;
t368 = t393 * t85;
t308 = t393 * t368;
t316 = pkin(5) * t368;
t326 = t185 * t209;
t346 = t393 * t192;
t123 = 0.1e1 / t275 ^ 2;
t377 = t123 * t79;
t378 = t393 * t79;
t278 = t192 * t316;
t67 = (t278 - t378) * t176;
t25 = t189 * t296 + t190 * t258 - t191 * t257 + (-(t185 * t67 + (pkin(2) * (t183 * t199 * t368 + t326 * t378) * t211 - t183 * (t79 * t346 - t316) * t323) * t176) * t308 + (-t199 * t243 + (t183 * t325 + (t176 - t198) * t185) * t368) * t377) * t177;
t353 = t25 * t192;
t390 = t177 - 0.2e1;
t13 = t198 * t353 - t390 * t249;
t398 = 0.2e1 * t13;
t200 = cos(qJ(3,2));
t201 = cos(qJ(2,2));
t195 = sin(qJ(2,2));
t320 = t195 * t200;
t152 = pkin(2) * t320 - t201 * pkin(5);
t194 = sin(qJ(3,2));
t386 = pkin(2) * t194;
t274 = t152 * t183 + t185 * t386;
t394 = 0.1e1 / t274;
t342 = t394 * t209;
t288 = t394 * t342;
t187 = legFrame(2,2);
t171 = sin(t187);
t174 = cos(t187);
t147 = -t171 * t205 + t174 * t206;
t114 = t147 * t184 - t337;
t332 = t185 * t195;
t335 = t183 * t200;
t86 = ((-t147 * t332 - t201 * t204) * t184 - (t201 * t147 - t195 * t327) * t182) * t194 - t114 * t335;
t268 = t86 * t288;
t117 = t182 * t147 + t161;
t329 = t185 * t201;
t383 = pkin(2) * t200;
t80 = -(t114 * t329 - t195 * t117) * t383 - (t114 * t332 + t201 * t117) * pkin(5);
t248 = t80 * t268;
t178 = 0.1e1 / t200;
t212 = t200 ^ 2;
t179 = 0.1e1 / t212;
t284 = t178 * t342;
t241 = t183 * t80 * t284;
t322 = t194 * t195;
t218 = t185 * t322 + t335;
t321 = t194 * t201;
t111 = t182 * t321 + t218 * t184;
t344 = t394 * t178;
t286 = t174 * t344;
t255 = t111 * t286;
t287 = t171 * t344;
t256 = t111 * t287;
t108 = t218 * t182 - t184 * t321;
t295 = t108 * t344;
t364 = t394 * t86;
t306 = t394 * t364;
t315 = pkin(5) * t364;
t343 = t394 * t194;
t125 = 0.1e1 / t274 ^ 2;
t373 = t125 * t80;
t374 = t394 * t80;
t277 = t194 * t315;
t68 = (t277 - t374) * t178;
t26 = t189 * t295 + t190 * t256 - t191 * t255 + (-(t185 * t68 + (pkin(2) * (t183 * t201 * t364 + t326 * t374) * t212 - t183 * (t80 * t343 - t315) * t320) * t178) * t306 - (t201 * t241 + (-t183 * t322 + (-t178 + t200) * t185) * t364) * t373) * t179;
t352 = t26 * t194;
t389 = t179 - 0.2e1;
t14 = t200 * t352 - t389 * t248;
t397 = 0.2e1 * t14;
t202 = cos(qJ(3,1));
t203 = cos(qJ(2,1));
t197 = sin(qJ(2,1));
t317 = t197 * t202;
t153 = pkin(2) * t317 - t203 * pkin(5);
t196 = sin(qJ(3,1));
t385 = pkin(2) * t196;
t273 = t153 * t183 + t185 * t385;
t395 = 0.1e1 / t273;
t339 = t395 * t209;
t283 = t395 * t339;
t188 = legFrame(1,2);
t172 = sin(t188);
t175 = cos(t188);
t148 = -t172 * t205 + t175 * t206;
t115 = t148 * t184 - t337;
t331 = t185 * t197;
t334 = t183 * t202;
t87 = ((-t148 * t331 - t203 * t204) * t184 - (t203 * t148 - t197 * t327) * t182) * t196 - t115 * t334;
t266 = t87 * t283;
t118 = t182 * t148 + t161;
t328 = t185 * t203;
t382 = pkin(2) * t202;
t81 = -(t115 * t328 - t197 * t118) * t382 - (t115 * t331 + t203 * t118) * pkin(5);
t247 = t81 * t266;
t180 = 0.1e1 / t202;
t213 = t202 ^ 2;
t181 = 0.1e1 / t213;
t279 = t180 * t339;
t239 = t183 * t81 * t279;
t319 = t196 * t197;
t217 = t185 * t319 + t334;
t318 = t196 * t203;
t112 = t182 * t318 + t217 * t184;
t341 = t395 * t180;
t281 = t175 * t341;
t253 = t112 * t281;
t282 = t172 * t341;
t254 = t112 * t282;
t109 = t217 * t182 - t184 * t318;
t294 = t109 * t341;
t360 = t395 * t87;
t304 = t395 * t360;
t314 = pkin(5) * t360;
t340 = t395 * t196;
t127 = 0.1e1 / t273 ^ 2;
t369 = t127 * t81;
t370 = t395 * t81;
t276 = t196 * t314;
t69 = (t276 - t370) * t180;
t27 = t189 * t294 + t190 * t254 - t191 * t253 + (-(t185 * t69 + (pkin(2) * (t183 * t203 * t360 + t326 * t370) * t213 - t183 * (t81 * t340 - t314) * t317) * t180) * t304 + (-t203 * t239 + (t183 * t319 + (t180 - t202) * t185) * t360) * t369) * t181;
t351 = t27 * t196;
t388 = t181 - 0.2e1;
t15 = t202 * t351 - t388 * t247;
t396 = 0.2e1 * t15;
t101 = pkin(5) * (t182 * t199 + t184 * t333) + (-t182 * t193 + t184 * t330) * t384;
t392 = t101 * t173;
t391 = t101 * t170;
t381 = t182 * g(3);
t169 = t184 * g(3);
t380 = t393 * t25;
t252 = g(1) * t173 - g(2) * t170;
t143 = t252 * t184 - t381;
t338 = t182 * t185;
t149 = t183 * g(1) + g(2) * t338;
t150 = g(1) * t338 - t183 * g(2);
t157 = t185 * t169;
t207 = pkin(5) ^ 2;
t208 = pkin(2) ^ 2;
t290 = t176 * t346;
t64 = -pkin(5) * t79 * t290 + (t176 * t207 + t198 * t208) * t368;
t154 = pkin(5) * t193 + t199 * t384;
t224 = -t151 * t185 + t183 * t387;
t94 = t184 * t154 + t224 * t182;
t88 = t170 * t275 + t94 * t173;
t89 = -t94 * t170 + t173 * t275;
t97 = -t182 * t154 + t224 * t184;
t55 = (t189 * t97 + t190 * t89 + t191 * t88) * t393 + (t64 * t308 - t67 * t377) * t176;
t216 = t149 * t170 - t150 * t173 - t55 * t183 - t157;
t40 = t193 * t143 - t216 * t199;
t379 = t393 * t40;
t376 = t394 * t26;
t251 = g(1) * t174 - g(2) * t171;
t144 = t251 * t184 - t381;
t285 = t178 * t343;
t65 = -pkin(5) * t80 * t285 + (t178 * t207 + t200 * t208) * t364;
t155 = pkin(5) * t195 + t201 * t383;
t223 = -t152 * t185 + t183 * t386;
t95 = t184 * t155 + t223 * t182;
t90 = t171 * t274 + t95 * t174;
t91 = -t95 * t171 + t174 * t274;
t98 = -t182 * t155 + t223 * t184;
t56 = (t189 * t98 + t190 * t91 + t191 * t90) * t394 + (t65 * t306 - t68 * t373) * t178;
t215 = t149 * t171 - t150 * t174 - t56 * t183 - t157;
t41 = t195 * t144 - t215 * t201;
t375 = t394 * t41;
t372 = t395 * t27;
t250 = g(1) * t175 - g(2) * t172;
t145 = t250 * t184 - t381;
t280 = t180 * t340;
t66 = -pkin(5) * t81 * t280 + (t180 * t207 + t202 * t208) * t360;
t156 = pkin(5) * t197 + t203 * t382;
t222 = -t153 * t185 + t183 * t385;
t96 = t184 * t156 + t222 * t182;
t92 = t172 * t273 + t96 * t175;
t93 = -t96 * t172 + t175 * t273;
t99 = -t182 * t156 + t222 * t184;
t57 = (t189 * t99 + t190 * t93 + t191 * t92) * t395 + (t66 * t304 - t69 * t369) * t180;
t214 = t149 * t172 - t150 * t175 - t57 * t183 - t157;
t42 = t197 * t145 - t214 * t203;
t371 = t395 * t42;
t367 = t393 * t88;
t366 = t393 * t89;
t365 = t393 * t97;
t363 = t394 * t90;
t362 = t394 * t91;
t361 = t394 * t98;
t359 = t395 * t92;
t358 = t395 * t93;
t357 = t395 * t99;
t356 = t199 * t25;
t355 = t201 * t26;
t354 = t203 * t27;
t350 = t85 ^ 2 * t123;
t349 = t86 ^ 2 * t125;
t348 = t87 ^ 2 * t127;
t102 = (-t182 * t195 + t184 * t329) * t383 + pkin(5) * (t182 * t201 + t184 * t332);
t313 = t102 * t376;
t103 = (-t182 * t197 + t184 * t328) * t382 + pkin(5) * (t182 * t203 + t184 * t331);
t312 = t103 * t372;
t311 = t110 * t379;
t310 = t111 * t375;
t309 = t112 * t371;
t210 = 0.1e1 / pkin(2) ^ 2;
t307 = t123 * t210 * t79 ^ 2;
t305 = t125 * t210 * t80 ^ 2;
t303 = t127 * t210 * t81 ^ 2;
t302 = t177 * t350;
t301 = t179 * t349;
t300 = t181 * t348;
t104 = (t182 * t330 + t184 * t193) * t384 + (t182 * t333 - t184 * t199) * pkin(5);
t299 = t104 * t347;
t105 = (t182 * t329 + t184 * t195) * t383 + (t182 * t332 - t184 * t201) * pkin(5);
t298 = t105 * t344;
t106 = (t182 * t328 + t184 * t197) * t382 + (t182 * t331 - t184 * t203) * pkin(5);
t297 = t106 * t341;
t272 = t40 * t296;
t271 = t41 * t295;
t270 = t42 * t294;
t267 = t26 * t285;
t265 = t27 * t280;
t264 = t101 * t291;
t263 = t101 * t292;
t262 = t102 * t287;
t261 = t102 * t286;
t260 = t103 * t282;
t259 = t103 * t281;
t246 = t102 * t267;
t245 = t103 * t265;
t244 = t302 * t346;
t242 = t301 * t343;
t240 = t300 * t340;
t238 = t40 * t258;
t237 = t41 * t256;
t236 = t42 * t254;
t235 = t40 * t257;
t234 = t41 * t255;
t233 = t42 * t253;
t232 = (t307 + t350) * t177 * t193 - t356;
t231 = (t305 + t349) * t179 * t195 - t355;
t230 = (t303 + t348) * t181 * t197 - t354;
t229 = 0.2e1 * t249;
t228 = 0.2e1 * t248;
t227 = 0.2e1 * t247;
t226 = t102 * t242;
t225 = t103 * t240;
t142 = t250 * t182 + t169;
t141 = t251 * t182 + t169;
t140 = t252 * t182 + t169;
t121 = t203 * t145;
t120 = t201 * t144;
t119 = t199 * t143;
t75 = t388 * t348;
t74 = t389 * t349;
t73 = t390 * t350;
t54 = -t172 * g(1) - t175 * g(2) + t57;
t53 = -t171 * g(1) - t174 * g(2) + t56;
t52 = -t170 * g(1) - t173 * g(2) + t55;
t51 = t183 * t142 - t54 * t185;
t50 = t183 * t141 - t53 * t185;
t49 = t183 * t140 - t52 * t185;
t48 = t121 + (-t142 * t185 - t183 * t54) * t197;
t47 = t120 + (-t141 * t185 - t183 * t53) * t195;
t46 = t119 + (-t140 * t185 - t183 * t52) * t193;
t45 = t214 * t197 + t121;
t44 = t215 * t195 + t120;
t43 = t216 * t193 + t119;
t39 = (-t185 * t66 * t266 - (-t196 * t153 * t239 + t185 * (-t180 * t276 + t202 * t370)) * t81 * t283) * t181 + (t106 * t189 + (t172 * t190 - t175 * t191) * t103) * t279;
t38 = (-t185 * t65 * t268 - (-t194 * t152 * t241 + t185 * (-t178 * t277 + t200 * t374)) * t80 * t288) * t179 + (t105 * t189 + (t171 * t190 - t174 * t191) * t102) * t284;
t37 = (-t185 * t64 * t269 - (-t192 * t151 * t243 + t185 * (-t176 * t278 + t198 * t378)) * t79 * t293) * t177 + (t104 * t189 + (t170 * t190 - t173 * t191) * t101) * t289;
t36 = t180 * t303 + t39 * t196;
t35 = t178 * t305 + t38 * t194;
t34 = t176 * t307 + t37 * t192;
t33 = -t181 * t196 * t303 + t39 * t202;
t32 = -t179 * t194 * t305 + t38 * t200;
t31 = -t177 * t192 * t307 + t37 * t198;
t30 = t181 * t203 * t227 + t197 * t39;
t29 = t179 * t201 * t228 + t195 * t38;
t28 = t177 * t199 * t229 + t193 * t37;
t24 = t195 * t26 + t201 * t301;
t23 = t195 * t301 - t355;
t22 = t197 * t27 + t203 * t300;
t21 = t193 * t25 + t199 * t302;
t20 = t197 * t300 - t354;
t19 = t193 * t302 - t356;
t18 = (t180 * t227 + t351) * t196;
t17 = (t178 * t228 + t352) * t194;
t16 = (t176 * t229 + t353) * t192;
t12 = t51 * t196 + t48 * t202;
t11 = t48 * t196 - t202 * t51;
t10 = t50 * t194 + t47 * t200;
t9 = t47 * t194 - t200 * t50;
t8 = t49 * t192 + t46 * t198;
t7 = t46 * t192 - t198 * t49;
t6 = (t230 * t196 - t202 * t30) * t183 - t185 * t36;
t5 = (t231 * t194 - t200 * t29) * t183 - t185 * t35;
t4 = (t232 * t192 - t198 * t28) * t183 - t185 * t34;
t3 = (-t196 * t30 - t230 * t202) * t183 + t185 * t33;
t2 = (-t194 * t29 - t231 * t200) * t183 + t185 * t32;
t1 = (-t192 * t28 - t232 * t198) * t183 + t185 * t31;
t58 = [t54 * t359 + t53 * t363 + t52 * t367, -t25 * t257 - t253 * t27 - t255 * t26, -t235 - t234 - t233 + (-t19 * t367 - t20 * t359 - t23 * t363) * t183, -t43 * t257 - t44 * t255 - t45 * t253 + (-t21 * t367 - t22 * t359 - t24 * t363) * t183, -t16 * t257 - t17 * t255 - t18 * t253 + (t174 * t226 + t175 * t225 + t244 * t392) * t209, (-t259 * t75 - t261 * t74 - t264 * t73) * t209 - 0.2e1 * t13 * t257 - 0.2e1 * t14 * t255 - 0.2e1 * t15 * t253, -t34 * t257 - t35 * t255 - t36 * t253 + (-t174 * t246 - t175 * t245 - t264 * t353) * t209, -t31 * t257 - t32 * t255 - t33 * t253 + (-t174 * t313 - t175 * t312 - t380 * t392) * t209, (-t259 * t39 - t261 * t38 - t264 * t37) * t209, -t173 * t311 - t174 * t310 - t175 * t309 + t1 * t367 + t2 * t363 + t3 * t359 + (-t11 * t259 - t261 * t9 - t264 * t7) * t209, t192 * t235 + t194 * t234 + t196 * t233 + t4 * t367 + t5 * t363 + t6 * t359 + (-t10 * t261 - t12 * t259 - t264 * t8) * t209, t191 - g(1); t54 * t358 + t53 * t362 + t52 * t366, t25 * t258 + t254 * t27 + t256 * t26, t238 + t237 + t236 + (-t19 * t366 - t20 * t358 - t23 * t362) * t183, t43 * t258 + t44 * t256 + t45 * t254 + (-t21 * t366 - t22 * t358 - t24 * t362) * t183, t16 * t258 + t17 * t256 + t18 * t254 + (-t171 * t226 - t172 * t225 - t244 * t391) * t209, (t260 * t75 + t262 * t74 + t263 * t73) * t209 + t258 * t398 + t256 * t397 + t254 * t396, t34 * t258 + t35 * t256 + t36 * t254 + (t171 * t246 + t172 * t245 + t263 * t353) * t209, t31 * t258 + t32 * t256 + t33 * t254 + (t171 * t313 + t172 * t312 + t380 * t391) * t209, (t260 * t39 + t262 * t38 + t263 * t37) * t209, t170 * t311 + t171 * t310 + t172 * t309 + t1 * t366 + t2 * t362 + t3 * t358 + (t11 * t260 + t262 * t9 + t263 * t7) * t209, -t192 * t238 - t194 * t237 - t196 * t236 + t4 * t366 + t5 * t362 + t6 * t358 + (t10 * t262 + t12 * t260 + t263 * t8) * t209, t190 - g(2); t54 * t357 + t53 * t361 + t52 * t365, t25 * t296 + t26 * t295 + t27 * t294, t272 + t271 + t270 + (-t19 * t365 - t20 * t357 - t23 * t361) * t183, t43 * t296 + t44 * t295 + t45 * t294 + (-t21 * t365 - t22 * t357 - t24 * t361) * t183, t16 * t296 + t17 * t295 + t18 * t294 + (-t104 * t244 - t105 * t242 - t106 * t240) * t209, (t297 * t75 + t298 * t74 + t299 * t73) * t209 + t296 * t398 + t295 * t397 + t294 * t396, t34 * t296 + t35 * t295 + t36 * t294 + (t104 * t25 * t290 + t105 * t267 + t106 * t265) * t209, t31 * t296 + t32 * t295 + t33 * t294 + (t104 * t380 + t105 * t376 + t106 * t372) * t209, (t297 * t39 + t298 * t38 + t299 * t37) * t209, t1 * t365 + t107 * t379 + t108 * t375 + t109 * t371 + t2 * t361 + t3 * t357 + (t11 * t297 + t298 * t9 + t299 * t7) * t209, -t192 * t272 - t194 * t271 - t196 * t270 + t4 * t365 + t5 * t361 + t6 * t357 + (t10 * t298 + t12 * t297 + t299 * t8) * t209, t189 - g(3);];
tauX_reg  = t58;
