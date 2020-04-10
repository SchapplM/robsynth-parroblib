% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRR1G2P2A0
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRR1G2P2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:59
% EndTime: 2020-03-09 21:25:03
% DurationCPUTime: 3.59s
% Computational Cost: add. (57762->320), mult. (43476->482), div. (2799->10), fcn. (22815->71), ass. (0->239)
t349 = 2 * xDP(3);
t245 = 0.1e1 / pkin(3);
t348 = t245 / 0.2e1;
t218 = pkin(7) + qJ(3,1);
t190 = sin(t218);
t231 = sin(qJ(3,1));
t129 = pkin(1) * t190 + t231 * pkin(2);
t101 = 0.1e1 / t129;
t223 = legFrame(1,2);
t279 = qJ(1,1) + pkin(7);
t170 = t223 + t279;
t164 = qJ(3,1) + t170;
t140 = cos(t164);
t171 = -t223 + t279;
t165 = qJ(3,1) + t171;
t141 = cos(t165);
t84 = -t141 + t140;
t301 = t101 * t84;
t347 = t301 / 0.2e1;
t134 = sin(t164);
t135 = sin(t165);
t81 = t134 + t135;
t302 = t101 * t81;
t346 = t302 / 0.2e1;
t222 = legFrame(2,2);
t278 = qJ(1,2) + pkin(7);
t168 = t222 + t278;
t162 = qJ(3,2) + t168;
t138 = cos(t162);
t169 = -t222 + t278;
t163 = qJ(3,2) + t169;
t139 = cos(t163);
t83 = -t139 + t138;
t217 = pkin(7) + qJ(3,2);
t189 = sin(t217);
t229 = sin(qJ(3,2));
t128 = pkin(1) * t189 + t229 * pkin(2);
t99 = 0.1e1 / t128;
t309 = t83 * t99;
t345 = t309 / 0.2e1;
t221 = legFrame(3,2);
t277 = qJ(1,3) + pkin(7);
t166 = t221 + t277;
t160 = qJ(3,3) + t166;
t136 = cos(t160);
t167 = -t221 + t277;
t161 = qJ(3,3) + t167;
t137 = cos(t161);
t82 = -t137 + t136;
t216 = pkin(7) + qJ(3,3);
t188 = sin(t216);
t227 = sin(qJ(3,3));
t127 = pkin(1) * t188 + t227 * pkin(2);
t97 = 0.1e1 / t127;
t310 = t82 * t97;
t344 = t310 / 0.2e1;
t132 = sin(t162);
t133 = sin(t163);
t80 = t132 + t133;
t311 = t80 * t99;
t343 = t311 / 0.2e1;
t130 = sin(t160);
t131 = sin(t161);
t79 = t130 + t131;
t312 = t79 * t97;
t342 = t312 / 0.2e1;
t220 = cos(pkin(7));
t341 = pkin(2) * t220;
t192 = cos(t217);
t224 = xDDP(3);
t235 = cos(qJ(3,2));
t242 = xDP(2);
t243 = xDP(1);
t202 = qJ(1,2) + t222;
t203 = qJ(1,2) - t222;
t53 = (sin(t202) + sin(t203)) * pkin(1) + (sin(t168) + sin(t169)) * pkin(2) + t80 * pkin(3);
t56 = -t83 * pkin(3) + (-cos(t168) + cos(t169)) * pkin(2) + (-cos(t202) + cos(t203)) * pkin(1);
t198 = qJ(1,2) + t217;
t179 = cos(t198);
t236 = cos(qJ(1,2));
t77 = -t236 * pkin(1) - pkin(2) * cos(t278) - pkin(3) * t179;
t322 = (t56 * t242 - t53 * t243 + t349 * t77) * t99 * t245;
t268 = t322 / 0.2e1;
t226 = xDDP(1);
t288 = t226 / 0.2e1;
t225 = xDDP(2);
t289 = t225 / 0.2e1;
t298 = t179 * t99;
t100 = 0.1e1 / t128 ^ 2;
t51 = t179 * t349 + t83 * t242 + t80 * t243;
t306 = t100 * t51;
t45 = t51 * t99 / 0.2e1;
t39 = t45 + t268;
t332 = pkin(3) * t39;
t23 = t99 * t268 * t332 + t224 * t298 + t288 * t311 + t289 * t309 - (-t332 + (-pkin(1) * t192 - pkin(2) * t235) * t45) * t306 / 0.2e1;
t191 = cos(t216);
t233 = cos(qJ(3,3));
t200 = qJ(1,3) + t221;
t201 = qJ(1,3) - t221;
t52 = (sin(t200) + sin(t201)) * pkin(1) + (sin(t166) + sin(t167)) * pkin(2) + t79 * pkin(3);
t55 = -t82 * pkin(3) + (-cos(t166) + cos(t167)) * pkin(2) + (-cos(t200) + cos(t201)) * pkin(1);
t197 = qJ(1,3) + t216;
t178 = cos(t197);
t234 = cos(qJ(1,3));
t76 = -t234 * pkin(1) - pkin(2) * cos(t277) - pkin(3) * t178;
t321 = (t55 * t242 - t52 * t243 + t349 * t76) * t97 * t245;
t267 = t321 / 0.2e1;
t299 = t178 * t97;
t49 = t178 * t349 + t82 * t242 + t79 * t243;
t98 = 0.1e1 / t127 ^ 2;
t319 = t49 * t98;
t43 = t49 * t97 / 0.2e1;
t37 = t43 + t267;
t333 = pkin(3) * t37;
t22 = t97 * t267 * t333 + t224 * t299 + t288 * t312 + t289 * t310 - (-t333 + (-pkin(1) * t191 - pkin(2) * t233) * t43) * t319 / 0.2e1;
t193 = cos(t218);
t237 = cos(qJ(3,1));
t204 = qJ(1,1) + t223;
t205 = qJ(1,1) - t223;
t54 = (sin(t204) + sin(t205)) * pkin(1) + (sin(t170) + sin(t171)) * pkin(2) + t81 * pkin(3);
t57 = -t84 * pkin(3) + (-cos(t170) + cos(t171)) * pkin(2) + (-cos(t204) + cos(t205)) * pkin(1);
t199 = qJ(1,1) + t218;
t180 = cos(t199);
t238 = cos(qJ(1,1));
t78 = -t238 * pkin(1) - pkin(2) * cos(t279) - pkin(3) * t180;
t320 = (t57 * t242 - t54 * t243 + t349 * t78) * t101 * t245;
t266 = t320 / 0.2e1;
t291 = t101 * t180;
t102 = 0.1e1 / t129 ^ 2;
t50 = t180 * t349 + t84 * t242 + t81 * t243;
t300 = t102 * t50;
t44 = t50 * t101 / 0.2e1;
t38 = t44 + t266;
t331 = t38 * pkin(3);
t24 = t101 * t266 * t331 + t224 * t291 + t288 * t302 + t289 * t301 - (-t331 + (-pkin(1) * t193 - t237 * pkin(2)) * t44) * t300 / 0.2e1;
t340 = -0.2e1 * pkin(2);
t339 = 0.2e1 * pkin(2);
t337 = -g(1) / 0.2e1;
t336 = g(1) / 0.2e1;
t335 = -g(2) / 0.2e1;
t334 = g(2) / 0.2e1;
t330 = t130 / 0.2e1;
t329 = t132 / 0.2e1;
t328 = t134 / 0.2e1;
t327 = -t137 / 0.2e1;
t326 = -t139 / 0.2e1;
t325 = -t141 / 0.2e1;
t324 = pkin(1) * sin(pkin(7));
t323 = pkin(1) / 0.2e1;
t318 = t52 * t97;
t317 = t53 * t99;
t316 = t55 * t97;
t315 = t56 * t99;
t314 = t76 * t97;
t313 = t77 * t99;
t308 = t225 - g(2);
t307 = t226 - g(1);
t305 = t101 * t54;
t304 = t101 * t57;
t303 = t101 * t78;
t297 = t52 * t226;
t296 = t53 * t226;
t295 = t54 * t226;
t294 = t55 * t225;
t293 = t56 * t225;
t292 = t57 * t225;
t290 = t224 * t245;
t287 = pkin(3) * t339;
t280 = 0.2e1 * pkin(1);
t215 = pkin(1) ^ 2 + pkin(2) ^ 2;
t244 = pkin(3) ^ 2;
t34 = t43 + t267 / 0.2e1;
t275 = (t34 * t233 * t287 + t215 * t43 + t37 * t244 + (pkin(3) * t191 * t34 + t43 * t341) * t280) * t319;
t273 = t49 ^ 2 * t98 / 0.4e1;
t35 = t44 + t266 / 0.2e1;
t272 = (t35 * t237 * t287 + t215 * t44 + t38 * t244 + (pkin(3) * t193 * t35 + t44 * t341) * t280) * t300;
t36 = t45 + t268 / 0.2e1;
t271 = (t36 * t235 * t287 + t215 * t45 + t39 * t244 + (pkin(3) * t192 * t36 + t45 * t341) * t280) * t306;
t265 = t100 * t51 ^ 2 / 0.4e1;
t264 = t102 * t50 ^ 2 / 0.4e1;
t228 = sin(qJ(1,3));
t206 = sin(t221);
t209 = cos(t221);
t94 = t209 * g(1) - t206 * g(2);
t70 = g(3) * t228 - t94 * t234;
t230 = sin(qJ(1,2));
t207 = sin(t222);
t210 = cos(t222);
t95 = t210 * g(1) - t207 * g(2);
t71 = g(3) * t230 - t95 * t236;
t232 = sin(qJ(1,1));
t208 = sin(t223);
t211 = cos(t223);
t96 = t211 * g(1) - t208 * g(2);
t72 = g(3) * t232 - t96 * t238;
t260 = t34 * t267;
t259 = t35 * t266;
t258 = t36 * t268;
t257 = t290 * t314;
t256 = t290 * t313;
t255 = t290 * t303;
t187 = t220 * pkin(1) + pkin(2);
t254 = t37 * (t233 * t187 - t227 * t324 + pkin(3)) / (t227 * t187 + t233 * t324) * t321;
t253 = t38 * (t237 * t187 - t231 * t324 + pkin(3)) / (t231 * t187 + t237 * t324) * t320;
t252 = t39 * (t235 * t187 - t229 * t324 + pkin(3)) / (t229 * t187 + t235 * t324) * t322;
t251 = g(1) * t327 + g(2) * t330 + t131 * t335 + t136 * t337 + g(3) * sin(t197);
t250 = g(1) * t326 + g(2) * t329 + t133 * t335 + t138 * t337 + g(3) * sin(t198);
t249 = g(1) * t325 + g(2) * t328 + t135 * t335 + t140 * t337 + g(3) * sin(t199);
t248 = g(1) * t330 + g(2) * t327 + g(3) * t178 + t131 * t336 + t136 * t334;
t247 = g(1) * t329 + g(2) * t326 + g(3) * t179 + t133 * t336 + t138 * t334;
t246 = g(1) * t328 + g(2) * t325 + g(3) * t180 + t135 * t336 + t140 * t334;
t75 = g(3) * t238 + t96 * t232;
t74 = g(3) * t236 + t95 * t230;
t73 = g(3) * t234 + t94 * t228;
t66 = t307 * t208 + t308 * t211;
t65 = t307 * t207 + t308 * t210;
t64 = t307 * t206 + t308 * t209;
t21 = t24 * pkin(1) + t72;
t20 = pkin(1) * t23 + t71;
t19 = pkin(1) * t22 + t70;
t18 = (t229 * t265 + t23 * t235) * pkin(2) + (t189 * t265 + t192 * t23) * pkin(1) + t250;
t17 = (t22 * t233 + t227 * t273) * pkin(2) + (t188 * t273 + t191 * t22) * pkin(1) + t251;
t16 = (-t231 * t24 + t237 * t264) * pkin(2) + (-t190 * t24 + t193 * t264) * pkin(1) + t246;
t15 = (-t229 * t23 + t235 * t265) * pkin(2) + (-t189 * t23 + t192 * t265) * pkin(1) + t247;
t14 = (-t22 * t227 + t233 * t273) * pkin(2) + (-t188 * t22 + t191 * t273) * pkin(1) + t248;
t13 = (t231 * t264 + t237 * t24) * pkin(2) + (t190 * t264 + t193 * t24) * pkin(1) + t249;
t12 = t255 - t253 / 0.2e1 + (-t272 + (t292 - t295) * t101) * t348 + t24;
t11 = t256 - t252 / 0.2e1 + (-t271 + (t293 - t296) * t99) * t348 + t23;
t10 = t257 - t254 / 0.2e1 + (-t275 + (t294 - t297) * t97) * t348 + t22;
t9 = t255 / 0.2e1 - t253 / 0.4e1 + (-t272 / 0.2e1 + (-t295 / 0.2e1 + t292 / 0.2e1) * t101) * t348 + t24;
t8 = t256 / 0.2e1 - t252 / 0.4e1 + (-t271 / 0.2e1 + (-t296 / 0.2e1 + t293 / 0.2e1) * t99) * t348 + t23;
t7 = t257 / 0.2e1 - t254 / 0.4e1 + (-t275 / 0.2e1 + (-t297 / 0.2e1 + t294 / 0.2e1) * t97) * t348 + t22;
t6 = (-t231 * t259 + t9 * t237) * t339 + (-t190 * t259 + t9 * t193) * t280 + t249;
t5 = (-t229 * t258 + t8 * t235) * t339 + (-t189 * t258 + t8 * t192) * t280 + t250;
t4 = (-t227 * t260 + t7 * t233) * t339 + (-t188 * t260 + t7 * t191) * t280 + t251;
t3 = (t231 * t9 + t237 * t259) * t340 + (-t9 * t190 - t193 * t259) * t280 + t246;
t2 = (t229 * t8 + t235 * t258) * t340 + (-t8 * t189 - t192 * t258) * t280 + t247;
t1 = (t227 * t7 + t233 * t260) * t340 + (-t7 * t188 - t191 * t260) * t280 + t248;
t25 = [t22 * t342 + t23 * t343 + t24 * t346, t70 * t342 + t71 * t343 + t72 * t346, t73 * t342 + t74 * t343 + t75 * t346, t206 * t64 + t207 * t65 + t208 * t66 + (t19 * t312 + t20 * t311 + t21 * t302) * t323, t10 * t342 + t12 * t346 + t11 * t343 + (-t10 * t318 - t11 * t317 - t12 * t305) * t348, t6 * t346 + t4 * t342 + t5 * t343 + (-t13 * t305 - t17 * t318 - t18 * t317) * t348, t1 * t342 + t3 * t346 + t2 * t343 + (-t14 * t318 - t15 * t317 - t16 * t305) * t348, t307; t22 * t344 + t23 * t345 + t24 * t347, t70 * t344 + t71 * t345 + t72 * t347, t73 * t344 + t74 * t345 + t75 * t347, t209 * t64 + t210 * t65 + t211 * t66 + (t19 * t310 + t20 * t309 + t21 * t301) * t323, t10 * t344 + t12 * t347 + t11 * t345 + (t10 * t316 + t11 * t315 + t12 * t304) * t348, t6 * t347 + t4 * t344 + t5 * t345 + (t13 * t304 + t17 * t316 + t18 * t315) * t348, t1 * t344 + t3 * t347 + t2 * t345 + (t14 * t316 + t15 * t315 + t16 * t304) * t348, t308; t22 * t299 + t23 * t298 + t24 * t291, t72 * t291 + t71 * t298 + t70 * t299, t75 * t291 + t74 * t298 + t73 * t299, (t19 * t299 + t20 * t298 + t21 * t291) * pkin(1), t10 * t299 + t12 * t291 + t11 * t298 + (t10 * t314 + t11 * t313 + t12 * t303) * t245, t6 * t291 + t4 * t299 + t5 * t298 + (t13 * t303 + t17 * t314 + t18 * t313) * t245, t1 * t299 + t3 * t291 + t2 * t298 + (t14 * t314 + t15 * t313 + t16 * t303) * t245, t224 - g(3);];
tauX_reg  = t25;
