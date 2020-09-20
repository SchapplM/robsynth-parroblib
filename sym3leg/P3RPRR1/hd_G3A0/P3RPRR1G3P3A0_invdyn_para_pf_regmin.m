% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRR1G3P3A0
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
% tauX_reg [3x8]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRR1G3P3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:56
% EndTime: 2020-03-09 21:27:00
% DurationCPUTime: 3.98s
% Computational Cost: add. (57762->326), mult. (43473->489), div. (2799->10), fcn. (22812->71), ass. (0->244)
t219 = legFrame(3,2);
t272 = qJ(1,3) + pkin(7);
t167 = t219 + t272;
t168 = -t219 + t272;
t198 = qJ(1,3) + t219;
t199 = qJ(1,3) - t219;
t161 = qJ(3,3) + t167;
t134 = sin(t161);
t162 = qJ(3,3) + t168;
t135 = sin(t162);
t278 = t134 - t135;
t52 = (-sin(t199) + sin(t198)) * pkin(1) + (-sin(t168) + sin(t167)) * pkin(2) + t278 * pkin(3);
t220 = legFrame(2,2);
t273 = qJ(1,2) + pkin(7);
t169 = t220 + t273;
t170 = -t220 + t273;
t200 = qJ(1,2) + t220;
t201 = qJ(1,2) - t220;
t163 = qJ(3,2) + t169;
t136 = sin(t163);
t164 = qJ(3,2) + t170;
t137 = sin(t164);
t277 = t136 - t137;
t53 = (-sin(t201) + sin(t200)) * pkin(1) + (-sin(t170) + sin(t169)) * pkin(2) + t277 * pkin(3);
t221 = legFrame(1,2);
t274 = qJ(1,1) + pkin(7);
t171 = t221 + t274;
t172 = -t221 + t274;
t202 = qJ(1,1) + t221;
t203 = qJ(1,1) - t221;
t165 = qJ(3,1) + t171;
t138 = sin(t165);
t166 = qJ(3,1) + t172;
t139 = sin(t166);
t276 = t138 - t139;
t54 = (-sin(t203) + sin(t202)) * pkin(1) + (-sin(t172) + sin(t171)) * pkin(2) + t276 * pkin(3);
t243 = 0.1e1 / pkin(3);
t358 = t243 / 0.2e1;
t144 = cos(t165);
t145 = cos(t166);
t81 = t145 + t144;
t216 = pkin(7) + qJ(3,1);
t188 = sin(t216);
t229 = sin(qJ(3,1));
t127 = pkin(1) * t188 + t229 * pkin(2);
t99 = 0.1e1 / t127;
t300 = t81 * t99;
t357 = t300 / 0.2e1;
t142 = cos(t163);
t143 = cos(t164);
t80 = t143 + t142;
t215 = pkin(7) + qJ(3,2);
t187 = sin(t215);
t227 = sin(qJ(3,2));
t126 = pkin(1) * t187 + t227 * pkin(2);
t97 = 0.1e1 / t126;
t301 = t80 * t97;
t356 = t301 / 0.2e1;
t140 = cos(t161);
t141 = cos(t162);
t79 = t141 + t140;
t214 = pkin(7) + qJ(3,3);
t186 = sin(t214);
t225 = sin(qJ(3,3));
t125 = pkin(1) * t186 + t225 * pkin(2);
t95 = 0.1e1 / t125;
t302 = t79 * t95;
t355 = t302 / 0.2e1;
t303 = t276 * t99;
t354 = -t303 / 0.2e1;
t304 = t277 * t97;
t353 = -t304 / 0.2e1;
t305 = t278 * t95;
t352 = -t305 / 0.2e1;
t218 = cos(pkin(7));
t351 = pkin(2) * t218;
t194 = cos(t216);
t235 = cos(qJ(3,1));
t240 = xDP(2);
t241 = xDP(1);
t239 = xDP(3);
t271 = 2 * t239;
t57 = -t81 * pkin(3) + (-cos(t171) - cos(t172)) * pkin(2) + (-cos(t202) - cos(t203)) * pkin(1);
t197 = qJ(1,1) + t216;
t175 = sin(t197);
t230 = sin(qJ(1,1));
t75 = t230 * pkin(1) + pkin(2) * sin(t274) + pkin(3) * t175;
t318 = (t54 * t240 + t57 * t241 + t75 * t271) * t99 * t243;
t258 = t318 / 0.2e1;
t224 = xDDP(1);
t280 = t224 / 0.2e1;
t223 = xDDP(2);
t281 = t223 / 0.2e1;
t100 = 0.1e1 / t127 ^ 2;
t339 = -2 * t239;
t50 = t175 * t339 - t276 * t240 + t81 * t241;
t297 = t100 * t50;
t44 = t50 * t99 / 0.2e1;
t38 = t44 + t258;
t334 = t38 * pkin(3);
t347 = t99 * t258 * t334 + t280 * t300 - t281 * t303 - (-t334 + (-pkin(1) * t194 - t235 * pkin(2)) * t44) * t297 / 0.2e1;
t193 = cos(t215);
t233 = cos(qJ(3,2));
t56 = -t80 * pkin(3) + (-cos(t169) - cos(t170)) * pkin(2) + (-cos(t200) - cos(t201)) * pkin(1);
t196 = qJ(1,2) + t215;
t174 = sin(t196);
t228 = sin(qJ(1,2));
t74 = t228 * pkin(1) + pkin(2) * sin(t273) + pkin(3) * t174;
t317 = (t53 * t240 + t56 * t241 + t74 * t271) * t97 * t243;
t257 = t317 / 0.2e1;
t51 = t174 * t339 - t277 * t240 + t80 * t241;
t98 = 0.1e1 / t126 ^ 2;
t315 = t51 * t98;
t45 = t51 * t97 / 0.2e1;
t39 = t45 + t257;
t335 = pkin(3) * t39;
t346 = t97 * t257 * t335 + t280 * t301 - t281 * t304 - (-t335 + (-pkin(1) * t193 - pkin(2) * t233) * t45) * t315 / 0.2e1;
t192 = cos(t214);
t231 = cos(qJ(3,3));
t55 = -t79 * pkin(3) + (-cos(t167) - cos(t168)) * pkin(2) + (-cos(t198) - cos(t199)) * pkin(1);
t195 = qJ(1,3) + t214;
t173 = sin(t195);
t226 = sin(qJ(1,3));
t73 = t226 * pkin(1) + pkin(2) * sin(t272) + pkin(3) * t173;
t319 = (t52 * t240 + t55 * t241 + t73 * t271) * t95 * t243;
t259 = t319 / 0.2e1;
t49 = t173 * t339 - t278 * t240 + t79 * t241;
t96 = 0.1e1 / t125 ^ 2;
t316 = t49 * t96;
t43 = t49 * t95 / 0.2e1;
t37 = t43 + t259;
t336 = pkin(3) * t37;
t345 = t95 * t259 * t336 + t280 * t302 - t281 * t305 - (-t336 + (-pkin(1) * t192 - pkin(2) * t231) * t43) * t316 / 0.2e1;
t341 = -0.2e1 * pkin(2);
t340 = 0.2e1 * pkin(2);
t338 = g(1) / 0.2e1;
t337 = -g(2) / 0.2e1;
t333 = t135 / 0.2e1;
t332 = t137 / 0.2e1;
t331 = t139 / 0.2e1;
t330 = t140 / 0.2e1;
t329 = t142 / 0.2e1;
t328 = t144 / 0.2e1;
t321 = pkin(1) * sin(pkin(7));
t320 = pkin(1) / 0.2e1;
t314 = t52 * t95;
t313 = t53 * t97;
t312 = t54 * t99;
t311 = t55 * t95;
t310 = t56 * t97;
t309 = t57 * t99;
t308 = t73 * t95;
t307 = t74 * t97;
t306 = t75 * t99;
t299 = t223 - g(2);
t298 = t224 - g(1);
t296 = t173 * t95;
t295 = t174 * t97;
t294 = t175 * t99;
t222 = xDDP(3);
t293 = t222 * t95;
t292 = t222 * t97;
t291 = t222 * t99;
t290 = t243 * t73;
t289 = t243 * t74;
t288 = t243 * t75;
t287 = t52 * t223;
t286 = t53 * t223;
t285 = t54 * t223;
t284 = t55 * t224;
t283 = t56 * t224;
t282 = t57 * t224;
t232 = cos(qJ(1,3));
t204 = sin(t219);
t207 = cos(t219);
t92 = t207 * g(1) - t204 * g(2);
t67 = g(3) * t232 + t92 * t226;
t279 = pkin(3) * t340;
t275 = 0.2e1 * pkin(1);
t213 = pkin(1) ^ 2 + pkin(2) ^ 2;
t242 = pkin(3) ^ 2;
t34 = t43 + t259 / 0.2e1;
t270 = (t34 * t231 * t279 + t213 * t43 + t37 * t242 + (pkin(3) * t192 * t34 + t43 * t351) * t275) * t316;
t36 = t45 + t257 / 0.2e1;
t269 = (t36 * t233 * t279 + t213 * t45 + t39 * t242 + (pkin(3) * t193 * t36 + t45 * t351) * t275) * t315;
t263 = t49 ^ 2 * t96 / 0.4e1;
t262 = t51 ^ 2 * t98 / 0.4e1;
t35 = t44 + t258 / 0.2e1;
t261 = (t35 * t235 * t279 + t213 * t44 + t38 * t242 + (pkin(3) * t194 * t35 + t44 * t351) * t275) * t297;
t256 = t100 * t50 ^ 2 / 0.4e1;
t234 = cos(qJ(1,2));
t205 = sin(t220);
t208 = cos(t220);
t93 = t208 * g(1) - t205 * g(2);
t69 = g(3) * t234 + t93 * t228;
t236 = cos(qJ(1,1));
t206 = sin(t221);
t209 = cos(t221);
t94 = t209 * g(1) - t206 * g(2);
t71 = g(3) * t236 + t94 * t230;
t255 = t34 * t259;
t254 = t35 * t258;
t253 = t36 * t257;
t185 = t218 * pkin(1) + pkin(2);
t252 = (t231 * t185 - t225 * t321 + pkin(3)) * t37 / (t225 * t185 + t231 * t321) * t319;
t251 = (t233 * t185 - t227 * t321 + pkin(3)) * t39 / (t227 * t185 + t233 * t321) * t317;
t250 = (t235 * t185 - t229 * t321 + pkin(3)) * t38 / (t229 * t185 + t235 * t321) * t318;
t249 = g(1) * t333 + g(2) * t330 + t134 * t338 + t141 * t337 + g(3) * cos(t195);
t248 = g(1) * t332 + g(2) * t329 + t136 * t338 + t143 * t337 + g(3) * cos(t196);
t247 = g(1) * t331 + g(2) * t328 + t138 * t338 + t145 * t337 + g(3) * cos(t197);
t246 = g(1) * t330 + g(2) * t333 - g(3) * t173 + t134 * t337 + t141 * t338;
t245 = g(1) * t329 + g(2) * t332 - g(3) * t174 + t136 * t337 + t143 * t338;
t244 = g(1) * t328 + g(2) * t331 - g(3) * t175 + t138 * t337 + t145 * t338;
t72 = -g(3) * t230 + t94 * t236;
t70 = -g(3) * t228 + t93 * t234;
t68 = -g(3) * t226 + t92 * t232;
t66 = t298 * t206 + t299 * t209;
t65 = t298 * t205 + t299 * t208;
t64 = t298 * t204 + t299 * t207;
t24 = -t175 * t291 + t347;
t23 = -t174 * t292 + t346;
t22 = -t173 * t293 + t345;
t21 = t22 * pkin(1) + t67;
t20 = t24 * pkin(1) + t71;
t19 = pkin(1) * t23 + t69;
t18 = (-t229 * t24 + t235 * t256) * pkin(2) + (-t188 * t24 + t194 * t256) * pkin(1) + t244;
t17 = (-t227 * t23 + t233 * t262) * pkin(2) + (-t187 * t23 + t193 * t262) * pkin(1) + t245;
t16 = (-t22 * t225 + t231 * t263) * pkin(2) + (-t186 * t22 + t192 * t263) * pkin(1) + t246;
t15 = (t229 * t256 + t235 * t24) * pkin(2) + (t188 * t256 + t194 * t24) * pkin(1) + t247;
t14 = (t227 * t262 + t23 * t233) * pkin(2) + (t187 * t262 + t193 * t23) * pkin(1) + t248;
t13 = (t22 * t231 + t225 * t263) * pkin(2) + (t186 * t263 + t192 * t22) * pkin(1) + t249;
t12 = (-t175 + t288) * t291 - t250 / 0.2e1 + (-t261 + (t282 + t285) * t99) * t358 + t347;
t11 = (-t174 + t289) * t292 - t251 / 0.2e1 + (-t269 + (t283 + t286) * t97) * t358 + t346;
t10 = (-t173 + t290) * t293 - t252 / 0.2e1 + (-t270 + (t284 + t287) * t95) * t358 + t345;
t9 = (-t175 + t288 / 0.2e1) * t291 - t250 / 0.4e1 + (-t261 / 0.2e1 + (t282 / 0.2e1 + t285 / 0.2e1) * t99) * t358 + t347;
t8 = (-t174 + t289 / 0.2e1) * t292 - t251 / 0.4e1 + (-t269 / 0.2e1 + (t283 / 0.2e1 + t286 / 0.2e1) * t97) * t358 + t346;
t7 = (-t173 + t290 / 0.2e1) * t293 - t252 / 0.4e1 + (-t270 / 0.2e1 + (t284 / 0.2e1 + t287 / 0.2e1) * t95) * t358 + t345;
t6 = (t229 * t9 + t235 * t254) * t341 + (-t9 * t188 - t194 * t254) * t275 + t244;
t5 = (t227 * t8 + t233 * t253) * t341 + (-t8 * t187 - t193 * t253) * t275 + t245;
t4 = (t225 * t7 + t231 * t255) * t341 + (-t7 * t186 - t192 * t255) * t275 + t246;
t3 = (-t229 * t254 + t9 * t235) * t340 + (-t188 * t254 + t9 * t194) * t275 + t247;
t2 = (-t227 * t253 + t8 * t233) * t340 + (-t187 * t253 + t8 * t193) * t275 + t248;
t1 = (-t225 * t255 + t7 * t231) * t340 + (-t186 * t255 + t7 * t192) * t275 + t249;
t25 = [t22 * t355 + t23 * t356 + t24 * t357, t67 * t355 + t69 * t356 + t71 * t357, t68 * t355 + t70 * t356 + t72 * t357, t204 * t64 + t205 * t65 + t206 * t66 + (t19 * t301 + t20 * t300 + t21 * t302) * t320, t10 * t355 + t11 * t356 + t12 * t357 + (t10 * t311 + t11 * t310 + t12 * t309) * t358, t1 * t355 + t2 * t356 + t3 * t357 + (t13 * t311 + t14 * t310 + t15 * t309) * t358, t4 * t355 + t5 * t356 + t6 * t357 + (t16 * t311 + t17 * t310 + t18 * t309) * t358, t298; t22 * t352 + t23 * t353 + t24 * t354, t67 * t352 + t69 * t353 + t71 * t354, t68 * t352 + t70 * t353 + t72 * t354, t207 * t64 + t208 * t65 + t209 * t66 + (-t19 * t304 - t20 * t303 - t21 * t305) * t320, t10 * t352 + t11 * t353 + t12 * t354 + (t10 * t314 + t11 * t313 + t12 * t312) * t358, t1 * t352 + t2 * t353 + t3 * t354 + (t13 * t314 + t14 * t313 + t15 * t312) * t358, t4 * t352 + t5 * t353 + t6 * t354 + (t16 * t314 + t17 * t313 + t18 * t312) * t358, t299; -t22 * t296 - t23 * t295 - t24 * t294, -t71 * t294 - t69 * t295 - t67 * t296, -t72 * t294 - t70 * t295 - t68 * t296, (-t19 * t295 - t20 * t294 - t21 * t296) * pkin(1), -t10 * t296 - t11 * t295 - t12 * t294 + (t10 * t308 + t11 * t307 + t12 * t306) * t243, -t1 * t296 - t2 * t295 - t3 * t294 + (t13 * t308 + t14 * t307 + t15 * t306) * t243, -t4 * t296 - t5 * t295 - t6 * t294 + (t16 * t308 + t17 * t307 + t18 * t306) * t243, t222 - g(3);];
tauX_reg  = t25;
