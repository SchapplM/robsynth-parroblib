% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR1G2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [4x15]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P4PRRRR1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:53
% EndTime: 2020-08-07 10:58:59
% DurationCPUTime: 6.06s
% Computational Cost: add. (13138->425), mult. (29236->842), div. (6272->22), fcn. (22400->26), ass. (0->356)
t191 = sin(qJ(2,1));
t159 = 0.1e1 / t191;
t196 = cos(qJ(3,1));
t169 = 0.1e1 / t196;
t336 = t159 * t169;
t197 = cos(qJ(2,1));
t202 = xP(4);
t143 = sin(t202);
t144 = cos(t202);
t206 = koppelP(1,2);
t210 = koppelP(1,1);
t118 = t143 * t210 + t144 * t206;
t122 = -t143 * t206 + t144 * t210;
t198 = xDP(4);
t173 = t198 ^ 2;
t182 = xDDP(4);
t184 = xDDP(2);
t100 = -t173 * t118 + t122 * t182 + t184;
t185 = xDDP(1);
t104 = -t118 * t182 - t173 * t122 + t185;
t181 = legFrame(1,2);
t138 = sin(t181);
t324 = t191 * t196;
t142 = cos(t181);
t190 = sin(qJ(3,1));
t344 = t142 * t190;
t109 = t138 * t324 - t344;
t112 = t138 * t190 + t142 * t324;
t183 = xDDP(3);
t233 = t196 ^ 2;
t170 = 0.1e1 / t233;
t171 = t169 * t170;
t211 = 0.1e1 / pkin(2);
t332 = t171 * t211;
t160 = 0.1e1 / t191 ^ 2;
t199 = xDP(3);
t325 = t190 * t197;
t200 = xDP(2);
t201 = xDP(1);
t93 = t118 * t142 + t138 * t122;
t92 = t138 * t200 - t201 * t142 + t93 * t198;
t84 = -t196 * t199 - t92 * t325;
t80 = t84 ^ 2;
t386 = t160 * t80;
t88 = t92 ^ 2;
t52 = -t138 * g(1) - t142 * g(2) + (t197 * t183 + (t100 * t112 + t104 * t109) * t169 + (t88 + t386) * t332) * t159;
t42 = g(3) * t191 + t52 * t197;
t393 = t336 * t42;
t189 = sin(qJ(2,2));
t156 = 0.1e1 / t189;
t194 = cos(qJ(3,2));
t165 = 0.1e1 / t194;
t338 = t156 * t165;
t195 = cos(qJ(2,2));
t205 = koppelP(2,2);
t209 = koppelP(2,1);
t117 = t143 * t209 + t144 * t205;
t121 = -t143 * t205 + t144 * t209;
t103 = -t117 * t182 - t173 * t121 + t185;
t180 = legFrame(2,2);
t137 = sin(t180);
t326 = t189 * t194;
t141 = cos(t180);
t188 = sin(qJ(3,2));
t346 = t141 * t188;
t108 = t137 * t326 - t346;
t111 = t137 * t188 + t141 * t326;
t230 = t194 ^ 2;
t166 = 0.1e1 / t230;
t167 = t165 * t166;
t333 = t167 * t211;
t157 = 0.1e1 / t189 ^ 2;
t327 = t188 * t195;
t96 = t117 * t141 + t137 * t121;
t91 = t137 * t200 - t201 * t141 + t96 * t198;
t83 = -t194 * t199 - t91 * t327;
t79 = t83 ^ 2;
t387 = t157 * t79;
t87 = t91 ^ 2;
t99 = -t173 * t117 + t121 * t182 + t184;
t51 = -t137 * g(1) - t141 * g(2) + (t195 * t183 + (t103 * t108 + t111 * t99) * t165 + (t87 + t387) * t333) * t156;
t41 = g(3) * t189 + t51 * t195;
t392 = t338 * t41;
t187 = sin(qJ(2,3));
t153 = 0.1e1 / t187;
t192 = cos(qJ(3,3));
t161 = 0.1e1 / t192;
t340 = t153 * t161;
t193 = cos(qJ(2,3));
t204 = koppelP(3,2);
t208 = koppelP(3,1);
t116 = t143 * t208 + t144 * t204;
t120 = -t143 * t204 + t144 * t208;
t102 = -t116 * t182 - t173 * t120 + t185;
t179 = legFrame(3,2);
t136 = sin(t179);
t328 = t187 * t192;
t140 = cos(t179);
t186 = sin(qJ(3,3));
t348 = t140 * t186;
t107 = t136 * t328 - t348;
t110 = t136 * t186 + t140 * t328;
t227 = t192 ^ 2;
t162 = 0.1e1 / t227;
t163 = t161 * t162;
t334 = t163 * t211;
t154 = 0.1e1 / t187 ^ 2;
t329 = t186 * t193;
t95 = t116 * t140 + t136 * t120;
t90 = t136 * t200 - t201 * t140 + t95 * t198;
t82 = -t192 * t199 - t90 * t329;
t78 = t82 ^ 2;
t388 = t154 * t78;
t86 = t90 ^ 2;
t98 = -t173 * t116 + t120 * t182 + t184;
t50 = -t136 * g(1) - t140 * g(2) + (t193 * t183 + (t102 * t107 + t110 * t98) * t161 + (t86 + t388) * t334) * t153;
t40 = g(3) * t187 + t50 * t193;
t391 = t340 * t40;
t175 = sin(qJ(2,4));
t146 = 0.1e1 / t175;
t176 = cos(qJ(3,4));
t148 = 0.1e1 / t176;
t343 = t146 * t148;
t177 = cos(qJ(2,4));
t203 = koppelP(4,2);
t207 = koppelP(4,1);
t115 = t143 * t207 + t144 * t203;
t119 = -t143 * t203 + t144 * t207;
t101 = -t115 * t182 - t173 * t119 + t185;
t178 = legFrame(4,2);
t135 = sin(t178);
t330 = t175 * t176;
t139 = cos(t178);
t174 = sin(qJ(3,4));
t350 = t139 * t174;
t105 = t135 * t330 - t350;
t106 = t135 * t174 + t139 * t330;
t218 = t176 ^ 2;
t149 = 0.1e1 / t218;
t150 = t148 * t149;
t341 = t150 * t211;
t147 = 0.1e1 / t175 ^ 2;
t331 = t174 * t177;
t94 = t115 * t139 + t135 * t119;
t89 = t135 * t200 - t201 * t139 + t94 * t198;
t81 = -t176 * t199 - t89 * t331;
t77 = t81 ^ 2;
t389 = t147 * t77;
t85 = t89 ^ 2;
t97 = -t173 * t115 + t119 * t182 + t184;
t49 = -t135 * g(1) - t139 * g(2) + (t177 * t183 + (t101 * t105 + t106 * t97) * t148 + (t85 + t389) * t341) * t146;
t37 = g(3) * t175 + t49 * t177;
t390 = t343 * t37;
t212 = 0.1e1 / pkin(2) ^ 2;
t385 = t184 - g(2);
t384 = t185 - g(1);
t383 = t146 * t89;
t382 = t148 * t94;
t381 = t153 * t90;
t380 = t156 * t91;
t379 = t159 * t92;
t378 = t161 * t95;
t377 = t165 * t96;
t376 = t169 * t93;
t271 = t101 * t139 - t135 * t97;
t302 = t146 * t331;
t287 = t149 * t302;
t342 = t146 * t177;
t25 = (t271 * t287 + (-t146 * t183 + (-(-t174 * t175 * t89 + t81 * t342) * t147 * t81 - (-t174 * t81 + t177 * t89) * t383) * t341) * t148) * t211;
t375 = t176 * t25;
t270 = t102 * t140 - t136 * t98;
t300 = t153 * t329;
t286 = t162 * t300;
t339 = t153 * t193;
t26 = (t270 * t286 + (-t153 * t183 + (-(-t186 * t187 * t90 + t82 * t339) * t154 * t82 - (-t186 * t82 + t193 * t90) * t381) * t334) * t161) * t211;
t374 = t192 * t26;
t269 = t103 * t141 - t137 * t99;
t298 = t156 * t327;
t285 = t166 * t298;
t337 = t156 * t195;
t27 = (t269 * t285 + (-t156 * t183 + (-(-t188 * t189 * t91 + t83 * t337) * t157 * t83 - (-t188 * t83 + t195 * t91) * t380) * t333) * t165) * t211;
t373 = t194 * t27;
t268 = t100 * t138 - t104 * t142;
t296 = t159 * t325;
t284 = t170 * t296;
t335 = t159 * t197;
t28 = (-t268 * t284 + (-t159 * t183 + (-(-t190 * t191 * t92 + t84 * t335) * t160 * t84 - (-t190 * t84 + t197 * t92) * t379) * t332) * t169) * t211;
t372 = t196 * t28;
t371 = t212 * t77;
t370 = t212 * t78;
t369 = t212 * t79;
t368 = t212 * t80;
t367 = t212 * t85;
t366 = t212 * t86;
t365 = t212 * t87;
t364 = t212 * t88;
t315 = t174 * t367;
t61 = -t271 * t211 * t148 + t150 * t315;
t363 = t61 * t174;
t362 = t61 * t176;
t314 = t186 * t366;
t62 = -t270 * t211 * t161 + t163 * t314;
t361 = t62 * t186;
t360 = t62 * t192;
t313 = t188 * t365;
t63 = -t269 * t211 * t165 + t167 * t313;
t359 = t63 * t188;
t358 = t63 * t194;
t312 = t190 * t364;
t64 = t268 * t211 * t169 + t171 * t312;
t357 = t64 * t190;
t356 = t64 * t196;
t355 = t135 * t148;
t354 = t136 * t161;
t353 = t137 * t165;
t352 = t138 * t169;
t351 = t139 * t148;
t349 = t140 * t161;
t347 = t141 * t165;
t345 = t142 * t169;
t151 = 0.1e1 / t218 ^ 2;
t323 = t151 * t389;
t322 = t151 * t371;
t164 = 0.1e1 / t227 ^ 2;
t321 = t164 * t388;
t168 = 0.1e1 / t230 ^ 2;
t320 = t168 * t387;
t172 = 0.1e1 / t233 ^ 2;
t319 = t172 * t386;
t318 = t164 * t370;
t317 = t168 * t369;
t316 = t172 * t368;
t311 = t105 * t343;
t310 = t106 * t343;
t309 = t107 * t340;
t308 = t108 * t338;
t307 = t109 * t336;
t306 = t110 * t340;
t305 = t111 * t338;
t304 = t112 * t336;
t303 = t149 * t342;
t301 = t162 * t339;
t299 = t166 * t337;
t297 = t170 * t335;
t295 = t212 * t81 * t383;
t294 = t212 * t82 * t381;
t293 = t212 * t83 * t380;
t292 = t212 * t84 * t379;
t291 = t174 * t323;
t290 = t186 * t321;
t289 = t188 * t320;
t288 = t190 * t319;
t283 = t94 * t287;
t282 = t95 * t286;
t281 = t96 * t285;
t280 = t93 * t284;
t279 = t135 * t287;
t278 = t136 * t286;
t277 = t137 * t285;
t276 = t138 * t284;
t275 = t139 * t287;
t274 = t140 * t286;
t273 = t141 * t285;
t272 = t142 * t284;
t9 = t174 * t375 + (0.2e1 * t148 - t150) * t295;
t267 = -0.2e1 * t9 * t287;
t266 = 0.2e1 * t170 * t292;
t265 = 0.2e1 * t149 * t295;
t264 = 0.2e1 * t162 * t294;
t263 = 0.2e1 * t166 * t293;
t123 = t139 * g(1) - t135 * g(2);
t38 = g(3) * t177 - t49 * t175;
t262 = t123 * t176 + t38 * t174 - t302 * t37;
t124 = t140 * g(1) - t136 * g(2);
t43 = g(3) * t193 - t50 * t187;
t261 = t124 * t192 + t43 * t186 - t300 * t40;
t125 = t141 * g(1) - t137 * g(2);
t45 = g(3) * t195 - t51 * t189;
t260 = t125 * t194 + t45 * t188 - t298 * t41;
t126 = t142 * g(1) - t138 * g(2);
t47 = g(3) * t197 - t52 * t191;
t259 = t126 * t196 + t47 * t190 - t296 * t42;
t10 = t186 * t374 + (0.2e1 * t161 - t163) * t294;
t258 = -0.2e1 * t10 * t286;
t11 = t188 * t373 + (0.2e1 * t165 - t167) * t293;
t257 = -0.2e1 * t11 * t285;
t12 = t190 * t372 + (0.2e1 * t169 - t171) * t292;
t256 = -0.2e1 * t12 * t284;
t53 = -t149 * t315 + t362;
t255 = t53 * t287 - t25;
t54 = -t162 * t314 + t360;
t254 = t54 * t286 - t26;
t55 = -t166 * t313 + t358;
t253 = t55 * t285 - t27;
t56 = -t170 * t312 + t356;
t252 = t56 * t284 - t28;
t251 = t262 * t148;
t250 = t261 * t161;
t249 = t260 * t165;
t248 = t259 * t169;
t57 = t148 * t367 + t363;
t247 = -t148 * t25 + t57 * t303;
t58 = t161 * t366 + t361;
t246 = -t161 * t26 + t58 * t301;
t59 = t165 * t365 + t359;
t245 = -t165 * t27 + t59 * t299;
t60 = t169 * t364 + t357;
t244 = -t169 * t28 + t60 * t297;
t145 = t174 ^ 2;
t243 = -t145 * t37 * t303 - t148 * (-t123 * t174 + t38 * t176);
t152 = t186 ^ 2;
t242 = -t152 * t40 * t301 - t161 * (-t124 * t186 + t43 * t192);
t155 = t188 ^ 2;
t241 = -t155 * t41 * t299 - t165 * (-t125 * t188 + t45 * t194);
t158 = t190 ^ 2;
t240 = -t158 * t42 * t297 - t169 * (-t126 * t190 + t47 * t196);
t239 = t186 * t246;
t238 = t190 * t244;
t237 = t247 * t174;
t236 = t245 * t188;
t213 = t211 * t212;
t114 = -t143 * t182 - t144 * t173;
t113 = -t143 * t173 + t144 * t182;
t76 = (-t109 * t118 + t112 * t122) * t336;
t75 = (-t108 * t117 + t111 * t121) * t338;
t74 = (-t107 * t116 + t110 * t120) * t340;
t73 = (-t105 * t115 + t106 * t119) * t343;
t72 = (t170 * t88 + t319) * t212;
t71 = (t166 * t87 + t320) * t212;
t70 = (t162 * t86 + t321) * t212;
t69 = (t149 * t85 + t323) * t212;
t68 = (-0.2e1 * t170 + t172) * t160 * t368;
t67 = (-0.2e1 * t166 + t168) * t157 * t369;
t66 = (-0.2e1 * t162 + t164) * t154 * t370;
t65 = (-0.2e1 * t149 + t151) * t147 * t371;
t24 = -t159 * t316 + t197 * t28;
t23 = -t156 * t317 + t195 * t27;
t22 = -t153 * t318 + t193 * t26;
t21 = -t160 * t197 * t316 - t28 * t191;
t20 = -t157 * t195 * t317 - t27 * t189;
t19 = -t154 * t193 * t318 - t26 * t187;
t18 = -t146 * t322 + t177 * t25;
t17 = -t147 * t177 * t322 - t25 * t175;
t16 = t158 * t28 + t190 * t266;
t15 = t155 * t27 + t188 * t263;
t14 = t152 * t26 + t186 * t264;
t13 = t145 * t25 + t174 * t265;
t8 = (t72 * t190 - t356) * t191 - t197 * (t190 * t28 + t266);
t7 = (t71 * t188 - t358) * t189 - t195 * (t188 * t27 + t263);
t6 = (t70 * t186 - t360) * t187 - t193 * (t186 * t26 + t264);
t5 = (-t72 * t196 - t357) * t191 + (-0.2e1 * t190 * t171 * t292 + t372) * t197;
t4 = (-t71 * t194 - t359) * t189 + (-0.2e1 * t188 * t167 * t293 + t373) * t195;
t3 = (-t70 * t192 - t361) * t187 + (-0.2e1 * t186 * t163 * t294 + t374) * t193;
t2 = (t69 * t174 - t362) * t175 - t177 * (t174 * t25 + t265);
t1 = (-t69 * t176 - t363) * t175 + (-0.2e1 * t174 * t150 * t295 + t375) * t177;
t29 = [t52 * t307 + t51 * t308 + t50 * t309 + t49 * t311, (t25 * t275 + t26 * t274 + t27 * t273 + t28 * t272) * t211, t18 * t311 + t22 * t309 + t23 * t308 + t24 * t307 + (t42 * t272 + t41 * t273 + t40 * t274 + t37 * t275) * t211, t17 * t311 + t19 * t309 + t20 * t308 + t21 * t307 + (t272 * t47 + t273 * t45 + t274 * t43 + t275 * t38) * t211, (t139 * t291 + t140 * t290 + t141 * t289 + t142 * t288) * t213 + (t13 * t275 + t14 * t274 + t15 * t273 + t16 * t272) * t211, (0.2e1 * t10 * t274 + 0.2e1 * t11 * t273 + 0.2e1 * t12 * t272 + 0.2e1 * t275 * t9 - t68 * t345 - t67 * t347 - t66 * t349 - t65 * t351) * t211, (t244 * t344 + t245 * t346 + t246 * t348 + t247 * t350) * t211, (t139 * t255 + t140 * t254 + t141 * t253 + t142 * t252) * t211, (-t345 * t64 - t347 * t63 - t349 * t62 - t351 * t61) * t211, t1 * t311 + t3 * t309 + t4 * t308 + t5 * t307 + (-t259 * t345 - t260 * t347 - t261 * t349 - t262 * t351) * t211, t2 * t311 + t6 * t309 + t7 * t308 + t8 * t307 + (t139 * t243 + t140 * t242 + t141 * t241 + t142 * t240) * t211, 0, t114, -t113, t384; t52 * t304 + t51 * t305 + t50 * t306 + t49 * t310, (-t25 * t279 - t26 * t278 - t27 * t277 - t28 * t276) * t211, t18 * t310 + t22 * t306 + t23 * t305 + t24 * t304 + (-t276 * t42 - t277 * t41 - t278 * t40 - t279 * t37) * t211, t17 * t310 + t19 * t306 + t20 * t305 + t21 * t304 + (-t276 * t47 - t277 * t45 - t278 * t43 - t279 * t38) * t211, (-t135 * t291 - t136 * t290 - t137 * t289 - t138 * t288) * t213 + (-t13 * t279 - t14 * t278 - t15 * t277 - t16 * t276) * t211, (t135 * t267 + t136 * t258 + t137 * t257 + t138 * t256 + t68 * t352 + t67 * t353 + t66 * t354 + t65 * t355) * t211, (-t135 * t237 - t136 * t239 - t137 * t236 - t138 * t238) * t211, (-t135 * t255 - t136 * t254 - t137 * t253 - t138 * t252) * t211, (t352 * t64 + t353 * t63 + t354 * t62 + t355 * t61) * t211, t1 * t310 + t3 * t306 + t4 * t305 + t5 * t304 + (t135 * t251 + t136 * t250 + t137 * t249 + t138 * t248) * t211, t2 * t310 + t6 * t306 + t7 * t305 + t8 * t304 + (-t135 * t243 - t136 * t242 - t137 * t241 - t138 * t240) * t211, 0, t113, t114, t385; t52 * t335 + t51 * t337 + t50 * t339 + t49 * t342, (-t25 * t343 - t26 * t340 - t27 * t338 - t28 * t336) * t211, t18 * t342 + t22 * t339 + t23 * t337 + t24 * t335 + (-t390 - t391 - t392 - t393) * t211, t17 * t342 + t19 * t339 + t20 * t337 + t21 * t335 + (-t336 * t47 - t338 * t45 - t340 * t43 - t343 * t38) * t211, (-t13 * t343 - t14 * t340 - t15 * t338 - t16 * t336) * t211, 0.2e1 * (-t10 * t340 - t11 * t338 - t12 * t336 - t343 * t9) * t211, (-t336 * t60 - t338 * t59 - t340 * t58 - t343 * t57) * t211, (-t336 * t56 - t338 * t55 - t340 * t54 - t343 * t53) * t211, 0, t1 * t342 + t3 * t339 + t4 * t337 + t5 * t335 + (-t146 * t37 - t153 * t40 - t156 * t41 - t159 * t42) * t211, t2 * t342 + t6 * t339 + t7 * t337 + t8 * t335 + (t174 * t390 + t186 * t391 + t188 * t392 + t190 * t393) * t211, 0, 0, 0, t183 - g(3); t73 * t49 + t74 * t50 + t75 * t51 + t76 * t52, (-t25 * t283 - t26 * t282 - t27 * t281 - t28 * t280) * t211, t73 * t18 + t74 * t22 + t75 * t23 + t76 * t24 + (-t280 * t42 - t281 * t41 - t282 * t40 - t283 * t37) * t211, t73 * t17 + t74 * t19 + t75 * t20 + t76 * t21 + (-t280 * t47 - t281 * t45 - t282 * t43 - t283 * t38) * t211, (-t288 * t93 - t289 * t96 - t290 * t95 - t291 * t94) * t213 + (-t13 * t283 - t14 * t282 - t15 * t281 - t16 * t280) * t211, (t256 * t93 + t257 * t96 + t258 * t95 + t267 * t94 + t376 * t68 + t67 * t377 + t66 * t378 + t65 * t382) * t211, (-t236 * t96 - t237 * t94 - t238 * t93 - t239 * t95) * t211, (-t252 * t93 - t253 * t96 - t254 * t95 - t255 * t94) * t211, (t376 * t64 + t63 * t377 + t62 * t378 + t61 * t382) * t211, t73 * t1 + t74 * t3 + t75 * t4 + t76 * t5 + (t248 * t93 + t249 * t96 + t250 * t95 + t251 * t94) * t211, t73 * t2 + t74 * t6 + t75 * t7 + t76 * t8 + (-t240 * t93 - t241 * t96 - t242 * t95 - t243 * t94) * t211, t182, -t384 * t143 + t385 * t144, -t385 * t143 - t384 * t144, 0;];
tauX_reg  = t29;
