% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V1G1A0
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR8V1G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:50:06
% EndTime: 2020-08-06 16:50:15
% DurationCPUTime: 8.22s
% Computational Cost: add. (33640->403), mult. (84317->832), div. (4888->18), fcn. (83977->22), ass. (0->354)
t200 = sin(pkin(3));
t212 = sin(qJ(2,2));
t338 = t200 * t212;
t218 = cos(qJ(2,2));
t217 = cos(qJ(3,2));
t378 = pkin(2) * t217;
t181 = pkin(5) * t212 + t218 * t378;
t199 = sin(pkin(6));
t201 = cos(pkin(6));
t314 = t212 * t217;
t178 = pkin(2) * t314 - t218 * pkin(5);
t202 = cos(pkin(3));
t211 = sin(qJ(3,2));
t380 = pkin(2) * t211;
t233 = -t178 * t202 + t200 * t380;
t120 = t181 * t201 + t199 * t233;
t123 = t199 * t181 - t201 * t233;
t204 = legFrame(2,3);
t188 = sin(t204);
t191 = cos(t204);
t102 = t120 * t191 - t188 * t123;
t105 = t120 * t188 + t123 * t191;
t207 = xDDP(2);
t208 = xDDP(1);
t195 = 0.1e1 / t217;
t257 = t178 * t200 + t202 * t380;
t391 = 0.1e1 / t257;
t221 = xDP(2);
t222 = xDP(1);
t169 = -t199 * t222 + t201 * t221;
t170 = t199 * t221 + t201 * t222;
t129 = t188 * t169 + t170 * t191;
t323 = t202 * t221;
t165 = t212 * t323 + t218 * t222;
t322 = t202 * t222;
t166 = t212 * t322 - t218 * t221;
t335 = t200 * t217;
t99 = ((-t199 * t165 - t166 * t201) * t191 - (t165 * t201 - t199 * t166) * t188) * t211 - t129 * t335;
t367 = t391 * t99;
t298 = t391 * t367;
t135 = 0.1e1 / t257 ^ 2;
t126 = t169 * t191 - t188 * t170;
t325 = t202 * t218;
t329 = t202 * t212;
t93 = -(t126 * t212 + t129 * t325) * t378 - pkin(5) * (-t218 * t126 + t129 * t329);
t371 = t135 * t93;
t223 = pkin(5) ^ 2;
t224 = pkin(2) ^ 2;
t340 = t195 * t211;
t272 = t391 * t340;
t75 = -pkin(5) * t93 * t272 + (t195 * t223 + t217 * t224) * t367;
t302 = pkin(5) * t367;
t258 = t211 * t302;
t372 = t391 * t93;
t78 = (-t258 + t372) * t195;
t376 = (t298 * t75 + t371 * t78) * t195;
t206 = xDDP(3);
t390 = t206 - g(3);
t385 = -(t102 * t208 + t105 * t207) * t391 - t376 - t390;
t395 = t385 * t338;
t214 = sin(qJ(2,1));
t311 = t214 * t200;
t220 = cos(qJ(2,1));
t219 = cos(qJ(3,1));
t377 = pkin(2) * t219;
t182 = pkin(5) * t214 + t220 * t377;
t310 = t214 * t219;
t179 = pkin(2) * t310 - t220 * pkin(5);
t213 = sin(qJ(3,1));
t381 = pkin(2) * t200;
t232 = -t179 * t202 + t213 * t381;
t121 = t182 * t201 + t199 * t232;
t124 = t199 * t182 - t201 * t232;
t205 = legFrame(1,3);
t189 = sin(t205);
t192 = cos(t205);
t103 = t121 * t192 - t189 * t124;
t106 = t121 * t189 + t124 * t192;
t328 = t202 * t213;
t150 = pkin(2) * t328 + t179 * t200;
t146 = 0.1e1 / t150;
t130 = t189 * t169 + t170 * t192;
t167 = t214 * t323 + t220 * t222;
t168 = t214 * t322 - t220 * t221;
t333 = t200 * t219;
t100 = ((-t199 * t167 - t168 * t201) * t192 - (t167 * t201 - t199 * t168) * t189) * t213 - t130 * t333;
t197 = 0.1e1 / t219;
t266 = t200 * t310;
t332 = t200 * t220;
t141 = -pkin(5) * t332 + (t266 + t328) * pkin(2);
t388 = 0.1e1 / t141;
t127 = t169 * t192 - t189 * t170;
t324 = t202 * t220;
t327 = t202 * t214;
t94 = -(t127 * t214 + t130 * t324) * t377 - pkin(5) * (-t220 * t127 + t130 * t327);
t369 = t388 * t94;
t296 = t388 * t369;
t342 = t146 * t197;
t339 = t197 * t213;
t269 = t388 * t339;
t360 = t100 * t146;
t76 = -pkin(5) * t94 * t269 + (t197 * t223 + t219 * t224) * t360;
t301 = pkin(5) * t360;
t256 = t213 * t301;
t79 = (t256 - t369) * t197;
t230 = t100 * t342 * t388 * t76 - t79 * t197 * t296;
t386 = -(t103 * t208 + t106 * t207) * t146 - t230 - t390;
t394 = t386 * t311;
t210 = sin(qJ(2,3));
t318 = t210 * t200;
t216 = cos(qJ(2,3));
t215 = cos(qJ(3,3));
t379 = pkin(2) * t215;
t180 = pkin(5) * t210 + t216 * t379;
t317 = t210 * t215;
t177 = pkin(2) * t317 - t216 * pkin(5);
t209 = sin(qJ(3,3));
t234 = -t177 * t202 + t209 * t381;
t119 = t180 * t201 + t199 * t234;
t122 = t199 * t180 - t201 * t234;
t203 = legFrame(3,3);
t187 = sin(t203);
t190 = cos(t203);
t101 = t119 * t190 - t187 * t122;
t104 = t119 * t187 + t122 * t190;
t331 = t202 * t209;
t148 = pkin(2) * t331 + t177 * t200;
t142 = 0.1e1 / t148;
t193 = 0.1e1 / t215;
t267 = t200 * t317;
t336 = t200 * t216;
t139 = -pkin(5) * t336 + (t267 + t331) * pkin(2);
t389 = 0.1e1 / t139;
t125 = t169 * t190 - t187 * t170;
t128 = t187 * t169 + t170 * t190;
t326 = t202 * t216;
t330 = t202 * t210;
t92 = -(t125 * t210 + t128 * t326) * t379 - pkin(5) * (-t216 * t125 + t128 * t330);
t374 = t389 * t92;
t300 = t389 * t374;
t344 = t142 * t193;
t341 = t193 * t209;
t275 = t389 * t341;
t163 = t210 * t323 + t216 * t222;
t164 = t210 * t322 - t216 * t221;
t337 = t200 * t215;
t98 = ((-t199 * t163 - t164 * t201) * t190 - (t163 * t201 - t199 * t164) * t187) * t209 - t128 * t337;
t368 = t142 * t98;
t74 = -pkin(5) * t92 * t275 + (t193 * t223 + t215 * t224) * t368;
t303 = pkin(5) * t368;
t259 = t209 * t303;
t77 = (t259 - t374) * t193;
t231 = t344 * t389 * t74 * t98 - t77 * t193 * t300;
t387 = -(t101 * t208 + t104 * t207) * t142 - t231 - t390;
t393 = t387 * t318;
t392 = t391 * t195;
t227 = t215 ^ 2;
t194 = 0.1e1 / t227;
t384 = t194 - 0.2e1;
t228 = t217 ^ 2;
t196 = 0.1e1 / t228;
t383 = t196 - 0.2e1;
t229 = t219 ^ 2;
t198 = 0.1e1 / t229;
t382 = t198 - 0.2e1;
t225 = 0.1e1 / pkin(2);
t351 = t389 * t225;
t274 = t193 * t351;
t245 = t200 * t92 * t274;
t154 = t201 * t187 + t190 * t199;
t158 = t199 * t216 + t201 * t330;
t161 = -t199 * t330 + t201 * t216;
t116 = (-t187 * t158 + t161 * t190) * t209 - t154 * t337;
t279 = t116 * t344;
t151 = -t199 * t187 + t190 * t201;
t113 = (-t158 * t190 - t187 * t161) * t209 - t151 * t337;
t282 = t113 * t344;
t321 = t202 * t225;
t352 = t389 * t209;
t33 = t207 * t279 + t208 * t282 + (-(t202 * t77 + (pkin(2) * (t321 * t374 + t336 * t368) * t227 - (t92 * t352 - t303) * t267) * t193) * t389 * t368 - (t216 * t245 + (-t209 * t318 + (-t193 + t215) * t202) * t368) * t300) * t194;
t375 = t389 * t33;
t348 = t391 * t225;
t271 = t195 * t348;
t243 = t200 * t93 * t271;
t155 = t201 * t188 + t191 * t199;
t159 = t199 * t218 + t201 * t329;
t162 = -t199 * t329 + t201 * t218;
t117 = (-t188 * t159 + t162 * t191) * t211 - t155 * t335;
t278 = t117 * t392;
t152 = -t199 * t188 + t191 * t201;
t114 = (-t159 * t191 - t188 * t162) * t211 - t152 * t335;
t281 = t114 * t392;
t334 = t200 * t218;
t349 = t391 * t211;
t36 = t207 * t278 + t208 * t281 + (-(-t202 * t78 + (pkin(2) * (t321 * t372 + t334 * t367) * t228 - t200 * (t93 * t349 - t302) * t314) * t195) * t298 - (t218 * t243 + (-t211 * t338 + (-t195 + t217) * t202) * t367) * t371) * t196;
t373 = t391 * t36;
t345 = t388 * t225;
t268 = t197 * t345;
t241 = t200 * t94 * t268;
t156 = t201 * t189 + t192 * t199;
t157 = t199 * t327 - t201 * t220;
t160 = t199 * t220 + t201 * t327;
t118 = (-t157 * t192 - t189 * t160) * t213 - t156 * t333;
t277 = t118 * t342;
t153 = -t199 * t189 + t192 * t201;
t115 = (t189 * t157 - t160 * t192) * t213 - t153 * t333;
t280 = t115 * t342;
t346 = t388 * t213;
t34 = t207 * t277 + t208 * t280 + (-(t202 * t79 + (pkin(2) * (t321 * t369 + t332 * t360) * t229 - (t94 * t346 - t301) * t266) * t197) * t388 * t360 - (t220 * t241 + (-t213 * t311 + (-t197 + t219) * t202) * t360) * t296) * t198;
t370 = t388 * t34;
t366 = t33 * t209;
t365 = t34 * t213;
t364 = t36 * t211;
t363 = t98 ^ 2 / t148 ^ 2;
t362 = t99 ^ 2 * t135;
t361 = t100 ^ 2 / t150 ^ 2;
t359 = t101 * t142;
t358 = t102 * t391;
t357 = t103 * t146;
t356 = t104 * t142;
t355 = t105 * t391;
t354 = t106 * t146;
t353 = t389 * t193;
t347 = t388 * t197;
t320 = t210 * t151;
t319 = t210 * t154;
t316 = t212 * t152;
t315 = t212 * t155;
t313 = t214 * t153;
t312 = t214 * t156;
t309 = t216 * t151;
t308 = t216 * t154;
t307 = t218 * t152;
t306 = t218 * t155;
t305 = t220 * t153;
t304 = t220 * t156;
t226 = 0.1e1 / pkin(2) ^ 2;
t299 = 0.1e1 / t139 ^ 2 * t226 * t92 ^ 2;
t297 = t135 * t226 * t93 ^ 2;
t295 = 0.1e1 / t141 ^ 2 * t226 * t94 ^ 2;
t183 = t199 * g(1) - t201 * g(2);
t184 = t201 * g(1) + t199 * g(2);
t240 = t183 * t190 + t184 * t187;
t49 = (-t200 * t387 + t202 * t240) * t216 + (-t187 * t183 + t184 * t190) * t210;
t294 = t49 * t341;
t239 = t183 * t191 + t184 * t188;
t57 = (-t200 * t385 + t202 * t239) * t218 + (-t188 * t183 + t184 * t191) * t212;
t293 = t57 * t340;
t238 = t183 * t192 + t184 * t189;
t50 = (-t200 * t386 + t202 * t238) * t220 + (-t189 * t183 + t184 * t192) * t214;
t292 = t50 * t339;
t291 = t194 * t363;
t290 = t196 * t362;
t289 = t198 * t361;
t107 = -(t202 * t309 - t319) * t379 - pkin(5) * (t202 * t320 + t308);
t288 = t107 * t353;
t108 = -(t202 * t308 + t320) * t379 - (t202 * t319 - t309) * pkin(5);
t287 = t108 * t353;
t109 = -(t202 * t307 - t315) * t378 - pkin(5) * (t202 * t316 + t306);
t286 = t109 * t392;
t110 = -(t202 * t306 + t316) * t378 - (t202 * t315 - t307) * pkin(5);
t285 = t110 * t392;
t111 = -(t202 * t305 - t312) * t377 - pkin(5) * (t202 * t313 + t304);
t284 = t111 * t347;
t112 = -(t202 * t304 + t313) * t377 - (t202 * t312 - t305) * pkin(5);
t283 = t112 * t347;
t276 = t142 * t351;
t273 = t391 * t348;
t270 = t146 * t345;
t31 = t216 * t33;
t265 = (t299 + t363) * t194 * t210 - t31;
t35 = t218 * t36;
t264 = (t297 + t362) * t196 * t212 - t35;
t32 = t220 * t34;
t263 = (t295 + t361) * t198 * t214 - t32;
t255 = t98 * t276;
t249 = t92 * t255;
t262 = 0.2e1 * (t215 * t366 - t384 * t249) * t344;
t253 = t99 * t273;
t248 = t93 * t253;
t261 = 0.2e1 * (t217 * t364 - t383 * t248) * t392;
t250 = t100 * t270;
t247 = t94 * t250;
t260 = 0.2e1 * (t219 * t365 - t382 * t247) * t342;
t254 = t33 * t275;
t252 = t36 * t272;
t251 = t34 * t269;
t246 = t291 * t352;
t244 = t290 * t349;
t242 = t289 * t346;
t237 = 0.2e1 * t249;
t236 = 0.2e1 * t248;
t235 = 0.2e1 * t247;
t176 = t192 * g(1) + t189 * g(2);
t175 = t191 * g(1) + t188 * g(2);
t174 = t190 * g(1) + t187 * g(2);
t173 = t189 * g(1) - t192 * g(2);
t172 = t188 * g(1) - t191 * g(2);
t171 = t187 * g(1) - t190 * g(2);
t85 = t382 * t361;
t84 = t383 * t362;
t83 = t384 * t363;
t63 = t200 * t239 + t202 * t385;
t62 = t200 * t238 + t202 * t386;
t61 = t200 * t240 + t202 * t387;
t60 = t175 * (t199 * t325 + t201 * t212) + t172 * (-t199 * t212 + t201 * t325) - t385 * t334;
t59 = -t172 * t159 + t175 * t162 + t395;
t58 = (-t183 * t329 + t218 * t184) * t191 + (-t218 * t183 - t184 * t329) * t188 + t395;
t56 = t176 * (t199 * t324 + t201 * t214) + t173 * (-t199 * t214 + t201 * t324) - t386 * t332;
t55 = t174 * (t199 * t326 + t201 * t210) + t171 * (-t199 * t210 + t201 * t326) - t387 * t336;
t54 = -t176 * t157 - t173 * t160 + t394;
t53 = -t171 * t158 + t174 * t161 + t393;
t52 = (-t183 * t327 + t220 * t184) * t192 + (-t220 * t183 - t184 * t327) * t189 + t394;
t51 = (-t183 * t330 + t216 * t184) * t190 + (-t216 * t183 - t184 * t330) * t187 + t393;
t48 = (t111 * t208 + t112 * t207) * t268 + (-t202 * t76 * t250 - (-t213 * t179 * t241 + t202 * (-t197 * t256 + t219 * t369)) * t94 * t270) * t198;
t47 = (t109 * t208 + t110 * t207) * t271 + (-t202 * t75 * t253 - (-t211 * t178 * t243 + t202 * (-t195 * t258 + t217 * t372)) * t93 * t273) * t196;
t46 = (t107 * t208 + t108 * t207) * t274 + (-t202 * t74 * t255 - (-t209 * t177 * t245 + t202 * (-t193 * t259 + t215 * t374)) * t92 * t276) * t194;
t45 = t197 * t295 + t48 * t213;
t44 = t195 * t297 + t47 * t211;
t43 = t193 * t299 + t46 * t209;
t42 = -t198 * t213 * t295 + t48 * t219;
t41 = -t196 * t211 * t297 + t47 * t217;
t40 = -t194 * t209 * t299 + t46 * t215;
t39 = t198 * t220 * t235 + t214 * t48;
t38 = t196 * t218 * t236 + t212 * t47;
t37 = t194 * t216 * t237 + t210 * t46;
t30 = t212 * t36 + t218 * t290;
t29 = -t212 * t290 + t35;
t28 = t214 * t34 + t220 * t289;
t27 = t210 * t33 + t216 * t291;
t26 = -t214 * t289 + t32;
t25 = -t210 * t291 + t31;
t24 = (t195 * t236 + t364) * t211;
t23 = (t197 * t235 + t365) * t213;
t22 = (t193 * t237 + t366) * t209;
t18 = t63 * t211 + t58 * t217;
t17 = t58 * t211 - t63 * t217;
t16 = t62 * t213 + t52 * t219;
t15 = t52 * t213 - t62 * t219;
t14 = t61 * t209 + t51 * t215;
t13 = t51 * t209 - t61 * t215;
t12 = (-t211 * t38 - t264 * t217) * t200;
t11 = (t264 * t211 - t217 * t38) * t200;
t10 = (-t213 * t39 - t263 * t219) * t200;
t9 = (t263 * t213 - t219 * t39) * t200;
t8 = (-t209 * t37 - t265 * t215) * t200;
t7 = (t265 * t209 - t215 * t37) * t200;
t6 = -t202 * t44 + t11;
t5 = t202 * t41 + t12;
t4 = -t202 * t45 + t9;
t3 = -t202 * t43 + t7;
t2 = t202 * t42 + t10;
t1 = t202 * t40 + t8;
t19 = [-t357 * t386 - t358 * t385 - t359 * t387, t280 * t34 + t281 * t36 + t282 * t33, t55 * t282 + t60 * t281 + t56 * t280 + (t25 * t359 + t26 * t357 + t29 * t358) * t200, t53 * t282 + t59 * t281 + t54 * t280 + (-t27 * t359 - t28 * t357 - t30 * t358) * t200, t22 * t282 + t24 * t281 + t23 * t280 + (-t107 * t246 - t109 * t244 - t111 * t242) * t225, (t284 * t85 + t286 * t84 + t288 * t83) * t225 + t113 * t262 + t114 * t261 + t115 * t260, t43 * t282 + t44 * t281 + t45 * t280 + (t107 * t254 + t109 * t252 + t111 * t251) * t225, t40 * t282 + t41 * t281 + t42 * t280 + (t107 * t375 + t109 * t373 + t111 * t370) * t225, (t284 * t48 + t286 * t47 + t288 * t46) * t225, (t103 * t2 + t115 * t50) * t146 + (t102 * t5 + t114 * t57) * t391 + (t1 * t101 + t113 * t49) * t142 + (t13 * t288 + t15 * t284 + t17 * t286) * t225, (t103 * t4 - t115 * t292) * t146 + (t102 * t6 - t114 * t293) * t391 + (t101 * t3 - t113 * t294) * t142 + (t14 * t288 + t16 * t284 + t18 * t286) * t225, t208 - g(1); -t354 * t386 - t355 * t385 - t356 * t387, t277 * t34 + t278 * t36 + t279 * t33, t55 * t279 + t60 * t278 + t56 * t277 + (t25 * t356 + t26 * t354 + t29 * t355) * t200, t53 * t279 + t59 * t278 + t54 * t277 + (-t27 * t356 - t28 * t354 - t30 * t355) * t200, t22 * t279 + t24 * t278 + t23 * t277 + (-t108 * t246 - t110 * t244 - t112 * t242) * t225, (t283 * t85 + t285 * t84 + t287 * t83) * t225 + t116 * t262 + t117 * t261 + t118 * t260, t43 * t279 + t44 * t278 + t45 * t277 + (t108 * t254 + t110 * t252 + t112 * t251) * t225, t40 * t279 + t41 * t278 + t42 * t277 + (t108 * t375 + t110 * t373 + t112 * t370) * t225, (t283 * t48 + t285 * t47 + t287 * t46) * t225, (t106 * t2 + t118 * t50) * t146 + (t105 * t5 + t117 * t57) * t391 + (t1 * t104 + t116 * t49) * t142 + (t13 * t287 + t15 * t283 + t17 * t285) * t225, (t106 * t4 - t118 * t292) * t146 + (t105 * t6 - t117 * t293) * t391 + (t104 * t3 - t116 * t294) * t142 + (t14 * t287 + t16 * t283 + t18 * t285) * t225, t207 - g(2); -(3 * g(3)) + (3 * t206) + (t357 + t358 + t359) * t208 + (t354 + t355 + t356) * t207 + t230 + t231 + t376, 0, (t25 + t26 + t29) * t200, (-t27 - t28 - t30) * t200, 0, 0, 0, 0, 0, t10 + t12 + t8 + (t40 + t41 + t42) * t202, t11 + t7 + t9 + (-t43 - t44 - t45) * t202, t390;];
tauX_reg  = t19;
