% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR8V2G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:44
% EndTime: 2020-08-06 17:35:53
% DurationCPUTime: 9.61s
% Computational Cost: add. (57419->446), mult. (110990->845), div. (4326->17), fcn. (106566->22), ass. (0->346)
t232 = sin(pkin(4));
t246 = sin(qJ(2,1));
t354 = t232 * t246;
t252 = cos(qJ(2,1));
t255 = pkin(7) + pkin(6);
t217 = t252 * t255;
t197 = pkin(2) * t246 - t217;
t251 = cos(qJ(3,1));
t234 = cos(pkin(4));
t245 = sin(qJ(3,1));
t342 = t234 * t245;
t230 = t251 ^ 2;
t396 = pkin(3) * t230;
t138 = 0.1e1 / ((pkin(3) * t342 + t197 * t232) * t251 + pkin(2) * t342 + t354 * t396);
t239 = xDDP(2);
t240 = xDDP(1);
t218 = pkin(2) ^ 2 + t255 ^ 2;
t256 = pkin(3) ^ 2;
t231 = sin(pkin(8));
t233 = cos(pkin(8));
t253 = xDP(2);
t254 = xDP(1);
t184 = t231 * t253 + t233 * t254;
t185 = -t231 * t254 + t233 * t253;
t237 = legFrame(1,3);
t221 = sin(t237);
t224 = cos(t237);
t144 = t184 * t224 + t221 * t185;
t214 = t251 * pkin(3) + pkin(2);
t183 = t246 * t214 - t217;
t326 = t246 * t255;
t363 = (t214 * t252 + t326) * t234;
t123 = -t144 * t363 - (-t221 * t184 + t185 * t224) * t183;
t194 = t214 * t342;
t349 = t232 * t251;
t162 = t183 * t349 + t194;
t158 = 0.1e1 / t162;
t366 = t158 * t245;
t295 = t123 * t366;
t410 = 0.2e1 * pkin(2);
t325 = pkin(3) * t410;
t341 = t234 * t246;
t117 = ((-t184 * t341 + t252 * t185) * t224 - (t252 * t184 + t185 * t341) * t221) * t245 - t144 * t349;
t327 = t246 * t251;
t165 = pkin(3) * t327 + t197;
t153 = t165 * t349 + t194;
t149 = 0.1e1 / t153;
t384 = t117 * t149;
t286 = (-t255 * t295 + (t230 * t256 + t251 * t325 + t218) * t384) * t138 * t384;
t257 = 0.1e1 / pkin(3);
t381 = t123 * t158;
t298 = t257 * t381;
t200 = pkin(2) * t252 + t326;
t355 = t232 * t245;
t262 = pkin(3) * t355 - t197 * t234;
t141 = t231 * t200 - t262 * t233;
t168 = -t231 * t221 + t224 * t233;
t174 = t231 * t341 - t233 * t252;
t177 = t231 * t252 + t233 * t341;
t316 = pkin(2) * t355;
t416 = t200 * t233 + t262 * t231;
t108 = (-t221 * t174 + t177 * t224) * t396 + (t141 * t224 + t416 * t221) * t251 - t168 * t316;
t390 = t108 * t138;
t171 = t233 * t221 + t224 * t231;
t105 = -(t174 * t224 + t221 * t177) * t396 + (-t141 * t221 + t416 * t224) * t251 + t171 * t316;
t393 = t105 * t138;
t304 = t255 * t384;
t280 = t245 * t304;
t93 = t280 - t381;
t289 = t251 * t286 + t239 * t390 + (pkin(2) * t298 - t93 * t251) * t138 * t381 + t240 * t393;
t238 = xDDP(3);
t420 = t238 - g(3);
t403 = -t289 - t420;
t423 = t403 * t354;
t244 = sin(qJ(2,2));
t356 = t232 * t244;
t250 = cos(qJ(2,2));
t216 = t250 * t255;
t196 = pkin(2) * t244 - t216;
t249 = cos(qJ(3,2));
t243 = sin(qJ(3,2));
t344 = t234 * t243;
t229 = t249 ^ 2;
t397 = pkin(3) * t229;
t137 = 0.1e1 / ((pkin(3) * t344 + t196 * t232) * t249 + pkin(2) * t344 + t356 * t397);
t236 = legFrame(2,3);
t220 = sin(t236);
t223 = cos(t236);
t143 = t184 * t223 + t220 * t185;
t213 = t249 * pkin(3) + pkin(2);
t182 = t244 * t213 - t216;
t329 = t244 * t255;
t364 = (t213 * t250 + t329) * t234;
t122 = -t143 * t364 - (-t220 * t184 + t185 * t223) * t182;
t193 = t213 * t344;
t351 = t232 * t249;
t161 = t182 * t351 + t193;
t156 = 0.1e1 / t161;
t367 = t156 * t243;
t296 = t122 * t367;
t343 = t234 * t244;
t116 = ((-t184 * t343 + t250 * t185) * t223 - (t250 * t184 + t185 * t343) * t220) * t243 - t143 * t351;
t330 = t244 * t249;
t164 = pkin(3) * t330 + t196;
t152 = t164 * t351 + t193;
t147 = 0.1e1 / t152;
t385 = t116 * t147;
t287 = (-t255 * t296 + (t229 * t256 + t249 * t325 + t218) * t385) * t137 * t385;
t382 = t122 * t156;
t299 = t257 * t382;
t199 = pkin(2) * t250 + t329;
t357 = t232 * t243;
t263 = pkin(3) * t357 - t196 * t234;
t140 = t231 * t199 - t263 * t233;
t167 = -t231 * t220 + t223 * t233;
t173 = t231 * t343 - t233 * t250;
t176 = t231 * t250 + t233 * t343;
t317 = pkin(2) * t357;
t415 = t199 * t233 + t263 * t231;
t107 = (-t220 * t173 + t176 * t223) * t397 + (t140 * t223 + t415 * t220) * t249 - t167 * t317;
t391 = t107 * t137;
t170 = t233 * t220 + t223 * t231;
t104 = -(t173 * t223 + t220 * t176) * t397 + (-t140 * t220 + t415 * t223) * t249 + t170 * t317;
t394 = t104 * t137;
t306 = t255 * t385;
t281 = t243 * t306;
t92 = t281 - t382;
t290 = t249 * t287 + t239 * t391 + (pkin(2) * t299 - t92 * t249) * t137 * t382 + t240 * t394;
t404 = -t290 - t420;
t422 = t404 * t356;
t242 = sin(qJ(2,3));
t358 = t232 * t242;
t248 = cos(qJ(2,3));
t215 = t248 * t255;
t195 = pkin(2) * t242 - t215;
t247 = cos(qJ(3,3));
t241 = sin(qJ(3,3));
t346 = t234 * t241;
t228 = t247 ^ 2;
t398 = pkin(3) * t228;
t136 = 0.1e1 / ((pkin(3) * t346 + t195 * t232) * t247 + pkin(2) * t346 + t358 * t398);
t235 = legFrame(3,3);
t219 = sin(t235);
t222 = cos(t235);
t142 = t184 * t222 + t219 * t185;
t212 = t247 * pkin(3) + pkin(2);
t181 = t242 * t212 - t215;
t332 = t242 * t255;
t365 = (t212 * t248 + t332) * t234;
t121 = -t142 * t365 - (-t219 * t184 + t185 * t222) * t181;
t192 = t212 * t346;
t353 = t232 * t247;
t160 = t181 * t353 + t192;
t154 = 0.1e1 / t160;
t368 = t154 * t241;
t297 = t121 * t368;
t345 = t234 * t242;
t115 = ((-t184 * t345 + t248 * t185) * t222 - (t248 * t184 + t185 * t345) * t219) * t241 - t142 * t353;
t333 = t242 * t247;
t163 = pkin(3) * t333 + t195;
t151 = t163 * t353 + t192;
t145 = 0.1e1 / t151;
t386 = t115 * t145;
t288 = t136 * (-t255 * t297 + (t228 * t256 + t247 * t325 + t218) * t386) * t386;
t383 = t121 * t154;
t300 = t257 * t383;
t198 = pkin(2) * t248 + t332;
t359 = t232 * t241;
t264 = pkin(3) * t359 - t195 * t234;
t139 = t231 * t198 - t264 * t233;
t166 = -t231 * t219 + t222 * t233;
t172 = t231 * t345 - t233 * t248;
t175 = t231 * t248 + t233 * t345;
t318 = pkin(2) * t359;
t414 = t198 * t233 + t264 * t231;
t106 = (-t219 * t172 + t175 * t222) * t398 + (t139 * t222 + t414 * t219) * t247 - t166 * t318;
t392 = t106 * t136;
t169 = t233 * t219 + t222 * t231;
t103 = -(t172 * t222 + t219 * t175) * t398 + (-t139 * t219 + t414 * t222) * t247 + t169 * t318;
t395 = t103 * t136;
t308 = t255 * t386;
t282 = t241 * t308;
t91 = t282 - t383;
t291 = t247 * t288 + t239 * t392 + (pkin(2) * t300 - t91 * t247) * t136 * t383 + t240 * t395;
t405 = -t291 - t420;
t421 = t405 * t358;
t201 = -t231 * g(1) + t233 * g(2);
t202 = t233 * g(1) + t231 * g(2);
t273 = t201 * t222 - t202 * t219;
t258 = 0.1e1 / pkin(3) ^ 2;
t303 = t121 ^ 2 / t160 ^ 2 * t258;
t124 = -t166 * t353 - (t166 * t345 + t248 * t169) * t241;
t127 = -t169 * t353 - (-t248 * t166 + t169 * t345) * t241;
t309 = t248 * t386;
t352 = t232 * t248;
t43 = (-((t232 * t309 + t234 * t300) * t398 + ((-t297 + t308) * t242 + pkin(2) * t309) * t353 + t91 * t234) * t386 + t124 * t240 + t127 * t239 - (t300 * t352 + (t228 * t234 - t333 * t359 - t234) * t386) * t383) * t136;
t413 = (t201 * t219 + t202 * t222) * t242 - pkin(6) * t303 - (t232 * t405 + t273 * t234) * t248 + t43 * t410;
t271 = t201 * t223 - t202 * t220;
t302 = t122 ^ 2 / t161 ^ 2 * t258;
t125 = -t167 * t351 - (t167 * t343 + t250 * t170) * t243;
t128 = -t170 * t351 - (-t250 * t167 + t170 * t343) * t243;
t307 = t250 * t385;
t350 = t232 * t250;
t44 = (-((t232 * t307 + t234 * t299) * t397 + ((-t296 + t306) * t244 + pkin(2) * t307) * t351 + t92 * t234) * t385 + t125 * t240 + t128 * t239 - (t299 * t350 + (t229 * t234 - t330 * t357 - t234) * t385) * t382) * t137;
t412 = (t201 * t220 + t202 * t223) * t244 - pkin(6) * t302 - (t232 * t404 + t271 * t234) * t250 + t44 * t410;
t269 = t201 * t224 - t202 * t221;
t301 = t123 ^ 2 / t162 ^ 2 * t258;
t126 = -t168 * t349 - (t168 * t341 + t252 * t171) * t245;
t129 = -t171 * t349 - (-t252 * t168 + t171 * t341) * t245;
t305 = t252 * t384;
t348 = t232 * t252;
t45 = (-((t232 * t305 + t234 * t298) * t396 + ((-t295 + t304) * t246 + pkin(2) * t305) * t349 + t93 * t234) * t384 + t126 * t240 + t129 * t239 - (t298 * t348 + (t230 * t234 - t327 * t355 - t234) * t384) * t381) * t138;
t411 = (t201 * t221 + t202 * t224) * t246 - pkin(6) * t301 - (t232 * t403 + t269 * t234) * t252 + t45 * t410;
t285 = t145 * t300;
t276 = t115 * t285;
t406 = pkin(6) / 0.2e1;
t335 = t240 * t257;
t336 = t239 * t257;
t337 = t234 * t257;
t347 = t232 * t257;
t133 = -t166 * t365 + t181 * t169;
t371 = t133 * t154;
t130 = -t181 * t166 - t169 * t365;
t374 = t130 * t154;
t399 = pkin(2) * t257;
t76 = t335 * t371 + t336 * t374 - t288 * t337 - (-t234 * t282 + (-t163 * t241 * t347 + (t247 * t399 + t228) * t234) * t383) * t285;
t409 = -0.2e1 * pkin(2) * t276 - 0.2e1 * t76 * t406;
t284 = t147 * t299;
t275 = t116 * t284;
t134 = -t167 * t364 + t182 * t170;
t370 = t134 * t156;
t131 = -t182 * t167 - t170 * t364;
t373 = t131 * t156;
t77 = t335 * t370 + t336 * t373 - t287 * t337 - (-t234 * t281 + (-t164 * t243 * t347 + (t249 * t399 + t229) * t234) * t382) * t284;
t408 = -0.2e1 * pkin(2) * t275 - 0.2e1 * t77 * t406;
t283 = t149 * t298;
t274 = t117 * t283;
t135 = -t168 * t363 + t183 * t171;
t369 = t135 * t158;
t132 = -t183 * t168 - t171 * t363;
t372 = t132 * t158;
t78 = t335 * t369 + t336 * t372 - t286 * t337 - (-t234 * t280 + (-t165 * t245 * t347 + (t251 * t399 + t230) * t234) * t381) * t283;
t407 = -0.2e1 * pkin(2) * t274 - 0.2e1 * t78 * t406;
t402 = 0.2e1 * t228 - 0.1e1;
t401 = 0.2e1 * t229 - 0.1e1;
t400 = 0.2e1 * t230 - 0.1e1;
t389 = t115 ^ 2 / t151 ^ 2;
t388 = t116 ^ 2 / t152 ^ 2;
t387 = t117 ^ 2 / t153 ^ 2;
t380 = t124 * t136;
t379 = t125 * t137;
t378 = t126 * t138;
t377 = t127 * t136;
t376 = t128 * t137;
t375 = t129 * t138;
t340 = t234 * t248;
t339 = t234 * t250;
t338 = t234 * t252;
t334 = t241 * t247;
t331 = t243 * t249;
t328 = t245 * t251;
t321 = 0.2e1 * t136 * (t402 * t276 + t43 * t334);
t320 = 0.2e1 * t137 * (t401 * t275 + t44 * t331);
t319 = 0.2e1 * t138 * (t400 * t274 + t45 * t328);
t315 = t43 * t368;
t314 = t154 * t247 * t43;
t313 = t44 * t367;
t312 = t156 * t249 * t44;
t311 = t45 * t366;
t310 = t158 * t251 * t45;
t40 = t248 * t43;
t294 = (t303 + t389) * t242 - t40;
t41 = t250 * t44;
t293 = (t302 + t388) * t244 - t41;
t42 = t252 * t45;
t292 = (t301 + t387) * t246 - t42;
t279 = t154 * t334 * t389;
t278 = t156 * t331 * t388;
t277 = t158 * t328 * t387;
t267 = 0.2e1 * t276;
t266 = 0.2e1 * t275;
t265 = 0.2e1 * t274;
t191 = t224 * g(1) + t221 * g(2);
t190 = t223 * g(1) + t220 * g(2);
t189 = t222 * g(1) + t219 * g(2);
t188 = t221 * g(1) - t224 * g(2);
t187 = t220 * g(1) - t223 * g(2);
t186 = t219 * g(1) - t222 * g(2);
t90 = t400 * t387;
t89 = t401 * t388;
t88 = t402 * t389;
t66 = t269 * t232 - t234 * t403;
t65 = t271 * t232 - t234 * t404;
t64 = t273 * t232 - t234 * t405;
t63 = t191 * (t231 * t338 + t233 * t246) + t188 * (-t231 * t246 + t233 * t338) - t403 * t348;
t62 = t190 * (t231 * t339 + t233 * t244) + t187 * (-t231 * t244 + t233 * t339) - t404 * t350;
t61 = t189 * (t231 * t340 + t233 * t242) + t186 * (-t231 * t242 + t233 * t340) - t405 * t352;
t60 = -t191 * t174 - t188 * t177 + t423;
t59 = -t190 * t173 - t187 * t176 + t422;
t58 = -t189 * t172 - t186 * t175 + t421;
t57 = -t245 * t301 + t78 * t251;
t56 = t78 * t245 + t251 * t301;
t55 = -t243 * t302 + t77 * t249;
t54 = t77 * t243 + t249 * t302;
t53 = -t241 * t303 + t76 * t247;
t52 = t76 * t241 + t247 * t303;
t48 = t246 * t78 + t252 * t265;
t47 = t244 * t77 + t250 * t266;
t46 = t242 * t76 + t248 * t267;
t39 = t244 * t44 + t250 * t388;
t38 = -t244 * t388 + t41;
t37 = t246 * t45 + t252 * t387;
t36 = t242 * t43 + t248 * t389;
t35 = -t246 * t387 + t42;
t34 = -t242 * t389 + t40;
t33 = (t45 * t245 + t251 * t265) * t245;
t32 = (t44 * t243 + t249 * t266) * t243;
t31 = (t43 * t241 + t247 * t267) * t241;
t27 = (t201 * t341 + t252 * t202) * t224 + (t252 * t201 - t202 * t341) * t221 + t423 + pkin(2) * t387 - t45 * pkin(6);
t26 = (t201 * t343 + t250 * t202) * t223 + (t250 * t201 - t202 * t343) * t220 + t422 + pkin(2) * t388 - t44 * pkin(6);
t25 = (t201 * t345 + t248 * t202) * t222 + (t248 * t201 - t202 * t345) * t219 + t421 + pkin(2) * t389 - t43 * pkin(6);
t24 = (-t245 * t48 - t292 * t251) * t232;
t23 = (t292 * t245 - t251 * t48) * t232;
t22 = (-t243 * t47 - t293 * t249) * t232;
t21 = (t293 * t243 - t249 * t47) * t232;
t20 = (-t241 * t46 - t294 * t247) * t232;
t19 = (t294 * t241 - t247 * t46) * t232;
t18 = -t66 * t245 + t27 * t251;
t17 = t27 * t245 + t66 * t251;
t16 = -t65 * t243 + t26 * t249;
t15 = t26 * t243 + t65 * t249;
t14 = -t64 * t241 + t25 * t247;
t13 = t25 * t241 + t64 * t247;
t12 = t245 * t407 + t411 * t251;
t11 = -t411 * t245 + t251 * t407;
t10 = t243 * t408 + t412 * t249;
t9 = -t412 * t243 + t249 * t408;
t8 = t241 * t409 + t413 * t247;
t7 = -t413 * t241 + t247 * t409;
t6 = t234 * t57 + t24;
t5 = -t234 * t56 + t23;
t4 = t234 * t55 + t22;
t3 = -t234 * t54 + t21;
t2 = t234 * t53 + t20;
t1 = -t234 * t52 + t19;
t28 = [-t393 * t403 - t394 * t404 - t395 * t405, t45 * t378 + t44 * t379 + t43 * t380, t61 * t380 + t62 * t379 + t63 * t378 + (t34 * t395 + t35 * t393 + t38 * t394) * t232, t58 * t380 + t59 * t379 + t60 * t378 + (-t36 * t395 - t37 * t393 - t39 * t394) * t232, t31 * t380 + t32 * t379 + t33 * t378 + (-t133 * t279 - t134 * t278 - t135 * t277) * t257, (-t90 * t369 - t89 * t370 - t88 * t371) * t257 + t124 * t321 + t125 * t320 + t126 * t319, t52 * t380 + t54 * t379 + t56 * t378 + (t133 * t315 + t134 * t313 + t135 * t311) * t257, t53 * t380 + t55 * t379 + t57 * t378 + (t133 * t314 + t134 * t312 + t135 * t310) * t257, (t78 * t369 + t77 * t370 + t76 * t371) * t257, (t105 * t6 + t12 * t126) * t138 + (t10 * t125 + t104 * t4) * t137 + (t103 * t2 + t124 * t8) * t136 + (t13 * t371 + t15 * t370 + t17 * t369) * t257, (t105 * t5 + t11 * t126) * t138 + (t104 * t3 + t125 * t9) * t137 + (t1 * t103 + t124 * t7) * t136 + (t14 * t371 + t16 * t370 + t18 * t369) * t257, t240 - g(1); -t390 * t403 - t391 * t404 - t392 * t405, t45 * t375 + t44 * t376 + t43 * t377, t61 * t377 + t62 * t376 + t63 * t375 + (t34 * t392 + t35 * t390 + t38 * t391) * t232, t58 * t377 + t59 * t376 + t60 * t375 + (-t36 * t392 - t37 * t390 - t39 * t391) * t232, t31 * t377 + t32 * t376 + t33 * t375 + (-t130 * t279 - t131 * t278 - t132 * t277) * t257, (-t90 * t372 - t89 * t373 - t88 * t374) * t257 + t127 * t321 + t128 * t320 + t129 * t319, t52 * t377 + t54 * t376 + t56 * t375 + (t130 * t315 + t131 * t313 + t132 * t311) * t257, t53 * t377 + t55 * t376 + t57 * t375 + (t130 * t314 + t131 * t312 + t132 * t310) * t257, (t78 * t372 + t77 * t373 + t76 * t374) * t257, (t108 * t6 + t12 * t129) * t138 + (t10 * t128 + t107 * t4) * t137 + (t106 * t2 + t127 * t8) * t136 + (t13 * t374 + t15 * t373 + t17 * t372) * t257, (t108 * t5 + t11 * t129) * t138 + (t107 * t3 + t128 * t9) * t137 + (t1 * t106 + t127 * t7) * t136 + (t14 * t374 + t16 * t373 + t18 * t372) * t257, t239 - g(2); (3 * t238) - (3 * g(3)) + t289 + t290 + t291, 0, (t34 + t35 + t38) * t232, (-t36 - t37 - t39) * t232, 0, 0, 0, 0, 0, t20 + t22 + t24 + (t53 + t55 + t57) * t234, t19 + t21 + t23 + (-t52 - t54 - t56) * t234, t420;];
tauX_reg  = t28;
