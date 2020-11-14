% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR12V1G2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x14]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR12V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:41
% EndTime: 2020-08-06 18:24:46
% DurationCPUTime: 4.81s
% Computational Cost: add. (15559->390), mult. (25461->726), div. (3228->14), fcn. (18972->18), ass. (0->302)
t177 = legFrame(3,2);
t150 = sin(t177);
t153 = cos(t177);
t127 = t153 * g(1) - t150 * g(2);
t184 = sin(qJ(1,3));
t190 = cos(qJ(1,3));
t106 = g(3) * t190 + t127 * t184;
t183 = sin(qJ(3,3));
t353 = t183 * pkin(3);
t143 = qJ(2,3) + t353;
t134 = 0.1e1 / t143 ^ 2;
t167 = 0.1e1 / t183;
t314 = t134 * t167;
t366 = (-pkin(5) - pkin(6));
t166 = pkin(1) - t366;
t195 = xDP(3);
t142 = t195 * t166;
t159 = pkin(3) * t195;
t196 = xDP(2);
t197 = xDP(1);
t239 = t150 * t196 - t153 * t197;
t189 = cos(qJ(3,3));
t361 = pkin(3) * t190;
t254 = (t189 + 0.1e1) * (t189 - 0.1e1) * t361;
t115 = t150 * t197 + t153 * t196;
t286 = t189 * t115;
t173 = t189 ^ 2;
t375 = -t173 + 0.1e1;
t52 = ((t239 * qJ(2,3) + t142) * t190 + (qJ(2,3) * t195 - t239 * t166) * t184 + pkin(3) * t286) * t183 - t239 * t254 + qJ(2,3) * t286 + t375 * t184 * t159;
t97 = -t239 * t184 + t190 * t195;
t384 = 0.2e1 * t97 * t52 * t314 - t106;
t178 = legFrame(2,2);
t151 = sin(t178);
t154 = cos(t178);
t128 = t154 * g(1) - t151 * g(2);
t186 = sin(qJ(1,2));
t192 = cos(qJ(1,2));
t107 = g(3) * t192 + t128 * t186;
t185 = sin(qJ(3,2));
t352 = t185 * pkin(3);
t144 = qJ(2,2) + t352;
t137 = 0.1e1 / t144 ^ 2;
t169 = 0.1e1 / t185;
t310 = t137 * t169;
t238 = t151 * t196 - t154 * t197;
t191 = cos(qJ(3,2));
t359 = pkin(3) * t192;
t253 = (t191 + 0.1e1) * (t191 - 0.1e1) * t359;
t117 = t151 * t197 + t154 * t196;
t284 = t191 * t117;
t174 = t191 ^ 2;
t374 = -t174 + 0.1e1;
t53 = ((t238 * qJ(2,2) + t142) * t192 + (qJ(2,2) * t195 - t238 * t166) * t186 + pkin(3) * t284) * t185 - t238 * t253 + qJ(2,2) * t284 + t374 * t186 * t159;
t98 = -t238 * t186 + t192 * t195;
t383 = 0.2e1 * t98 * t53 * t310 - t107;
t179 = legFrame(1,2);
t152 = sin(t179);
t155 = cos(t179);
t129 = t155 * g(1) - t152 * g(2);
t188 = sin(qJ(1,1));
t194 = cos(qJ(1,1));
t108 = g(3) * t194 + t129 * t188;
t187 = sin(qJ(3,1));
t351 = t187 * pkin(3);
t145 = qJ(2,1) + t351;
t140 = 0.1e1 / t145 ^ 2;
t171 = 0.1e1 / t187;
t306 = t140 * t171;
t237 = t152 * t196 - t155 * t197;
t193 = cos(qJ(3,1));
t357 = pkin(3) * t194;
t252 = (t193 + 0.1e1) * (t193 - 0.1e1) * t357;
t119 = t152 * t197 + t155 * t196;
t282 = t193 * t119;
t175 = t193 ^ 2;
t373 = -t175 + 0.1e1;
t54 = ((t237 * qJ(2,1) + t142) * t194 + (qJ(2,1) * t195 - t237 * t166) * t188 + pkin(3) * t282) * t187 - t237 * t252 + qJ(2,1) * t282 + t373 * t188 * t159;
t99 = -t237 * t188 + t194 * t195;
t382 = 0.2e1 * t99 * t54 * t306 - t108;
t168 = 0.1e1 / t183 ^ 2;
t207 = 0.1e1 / pkin(3) ^ 2;
t326 = t115 ^ 2 * t207;
t266 = t168 * t326;
t170 = 0.1e1 / t185 ^ 2;
t325 = t117 ^ 2 * t207;
t264 = t170 * t325;
t172 = 0.1e1 / t187 ^ 2;
t324 = t119 ^ 2 * t207;
t262 = t172 * t324;
t121 = t127 * t190;
t198 = pkin(1) + pkin(5);
t156 = g(3) * t184;
t109 = t184 * t143 + t166 * t190;
t133 = 0.1e1 / t143;
t135 = t133 * t134;
t180 = xDDP(3);
t181 = xDDP(2);
t182 = xDDP(1);
t200 = qJ(2,3) ^ 2;
t205 = pkin(3) ^ 2;
t233 = -(pkin(6) ^ 2) - t205 + ((-2 * pkin(6) - pkin(5)) * pkin(5)) + ((2 * t366 - pkin(1)) * pkin(1));
t293 = t166 * t167;
t269 = t52 * t293;
t206 = 0.1e1 / pkin(3);
t290 = t167 * t206;
t316 = t133 * t189;
t323 = t115 * t133;
t341 = t166 * t97;
t349 = t134 * t97;
t131 = qJ(2,3) * t190 - t166 * t184;
t287 = t189 * qJ(2,3);
t362 = pkin(3) * t189;
t83 = (-t131 * t153 + t150 * t362) * t183 + t153 * t254 + t150 * t287;
t86 = (t131 * t150 + t153 * t362) * t183 + t375 * t150 * t361 + t153 * t287;
t236 = -t286 * t293 * t349 + t97 * t135 * t269 - ((t115 * t167 + t316 * t341) * t183 + qJ(2,3) * t115 * t290) * t168 * t323 + ((t269 + (-0.2e1 * qJ(2,3) * t353 + t173 * t205 - t200 + t233) * t97) * t349 - t109 * t180 + (-t86 * t181 - t83 * t182) * t167) * t133;
t229 = t156 + t236;
t94 = t97 ^ 2;
t350 = t134 * t94;
t224 = qJ(2,3) * t350 + t229;
t37 = (0.2e1 * t134 * t286 - t135 * t52) * t97 * t167 + (t190 * t180 - (t167 * t52 - t341) * t349 + (-t150 * t181 + t153 * t182) * t184) * t133;
t381 = -t198 * t37 + t121 - t224;
t122 = t128 * t192;
t157 = g(3) * t186;
t110 = t186 * t144 + t166 * t192;
t136 = 0.1e1 / t144;
t138 = t136 * t137;
t201 = qJ(2,2) ^ 2;
t292 = t166 * t169;
t268 = t53 * t292;
t289 = t169 * t206;
t312 = t136 * t191;
t322 = t117 * t136;
t340 = t166 * t98;
t346 = t137 * t98;
t132 = qJ(2,2) * t192 - t166 * t186;
t285 = t191 * qJ(2,2);
t360 = pkin(3) * t191;
t84 = (-t132 * t154 + t151 * t360) * t185 + t154 * t253 + t151 * t285;
t87 = (t132 * t151 + t154 * t360) * t185 + t374 * t151 * t359 + t154 * t285;
t235 = -t284 * t292 * t346 + t98 * t138 * t268 - ((t117 * t169 + t312 * t340) * t185 + qJ(2,2) * t117 * t289) * t170 * t322 + ((t268 + (-0.2e1 * qJ(2,2) * t352 + t174 * t205 - t201 + t233) * t98) * t346 - t110 * t180 + (-t87 * t181 - t84 * t182) * t169) * t136;
t228 = t157 + t235;
t95 = t98 ^ 2;
t347 = t137 * t95;
t225 = qJ(2,2) * t347 + t228;
t38 = (0.2e1 * t137 * t284 - t138 * t53) * t98 * t169 + (t192 * t180 - (t169 * t53 - t340) * t346 + (-t151 * t181 + t154 * t182) * t186) * t136;
t380 = -t198 * t38 + t122 - t225;
t123 = t129 * t194;
t158 = g(3) * t188;
t111 = t188 * t145 + t166 * t194;
t139 = 0.1e1 / t145;
t141 = t139 * t140;
t202 = qJ(2,1) ^ 2;
t291 = t166 * t171;
t267 = t54 * t291;
t288 = t171 * t206;
t308 = t139 * t193;
t321 = t119 * t139;
t339 = t166 * t99;
t343 = t140 * t99;
t130 = t194 * qJ(2,1) - t166 * t188;
t283 = t193 * qJ(2,1);
t358 = pkin(3) * t193;
t82 = (-t130 * t155 + t152 * t358) * t187 + t155 * t252 + t152 * t283;
t85 = (t130 * t152 + t155 * t358) * t187 + t373 * t152 * t357 + t155 * t283;
t234 = -t282 * t291 * t343 + t99 * t141 * t267 - ((t119 * t171 + t308 * t339) * t187 + qJ(2,1) * t119 * t288) * t172 * t321 + ((t267 + (-0.2e1 * qJ(2,1) * t351 + t175 * t205 - t202 + t233) * t99) * t343 - t111 * t180 + (-t85 * t181 - t82 * t182) * t171) * t139;
t227 = t158 + t234;
t96 = t99 ^ 2;
t344 = t140 * t96;
t226 = qJ(2,1) * t344 + t227;
t39 = (0.2e1 * t140 * t282 - t141 * t54) * t99 * t171 + (t194 * t180 - (t171 * t54 - t339) * t343 + (-t152 * t181 + t155 * t182) * t188) * t139;
t379 = -t198 * t39 + t123 - t226;
t363 = pkin(1) * t39;
t378 = t363 - t123;
t364 = pkin(1) * t38;
t377 = t364 - t122;
t365 = pkin(1) * t37;
t376 = t365 - t121;
t372 = t189 * t266;
t371 = t191 * t264;
t370 = t193 * t262;
t369 = 0.2e1 * qJ(2,1);
t368 = 0.2e1 * qJ(2,2);
t367 = 0.2e1 * qJ(2,3);
t147 = 0.2e1 * t173 - 0.1e1;
t148 = 0.2e1 * t174 - 0.1e1;
t149 = 0.2e1 * t175 - 0.1e1;
t348 = t135 * t94;
t345 = t138 * t95;
t342 = t141 * t96;
t338 = t167 * t83;
t337 = t167 * t86;
t336 = t169 * t84;
t335 = t169 * t87;
t334 = t171 * t82;
t333 = t171 * t85;
t332 = t189 * t37;
t331 = t191 * t38;
t330 = t193 * t39;
t91 = -t167 * t372 + (-t150 * t182 - t153 * t181) * t290;
t329 = t91 * t183;
t92 = -t169 * t371 + (-t151 * t182 - t154 * t181) * t289;
t328 = t92 * t185;
t93 = -t171 * t370 + (-t152 * t182 - t155 * t181) * t288;
t327 = t93 * t187;
t315 = t133 * t190;
t311 = t136 * t192;
t307 = t139 * t194;
t305 = t150 * t167;
t304 = t150 * t184;
t303 = t151 * t169;
t302 = t151 * t186;
t301 = t152 * t171;
t300 = t152 * t188;
t299 = t153 * t167;
t298 = t153 * t184;
t297 = t154 * t169;
t296 = t154 * t186;
t295 = t155 * t171;
t294 = t155 * t188;
t251 = t206 * t97 * t323;
t281 = (0.2e1 * t251 + t332) * t316;
t280 = t189 * t350;
t279 = t167 * t348;
t250 = t206 * t98 * t322;
t278 = (0.2e1 * t250 + t331) * t312;
t277 = t191 * t347;
t276 = t169 * t345;
t243 = t99 * t206 * t321;
t275 = (0.2e1 * t243 + t330) * t308;
t274 = t193 * t344;
t273 = t171 * t342;
t272 = t167 * t332;
t271 = t169 * t331;
t270 = t171 * t330;
t260 = t133 * t304;
t259 = t133 * t298;
t258 = t136 * t302;
t257 = t136 * t296;
t256 = t139 * t300;
t255 = t139 * t294;
t249 = t184 * t281;
t248 = t147 * t94 * t314;
t247 = t186 * t278;
t246 = t148 * t95 * t310;
t245 = t188 * t275;
t244 = t149 * t96 * t306;
t73 = -t167 * t326 + t91 * t189;
t74 = -t169 * t325 + t92 * t191;
t75 = -t171 * t324 + t93 * t193;
t242 = t167 * t251;
t241 = t169 * t250;
t240 = t171 * t243;
t25 = t37 * t367 + t384;
t26 = t38 * t368 + t383;
t27 = t39 * t369 + t382;
t199 = pkin(1) * g(3);
t126 = t152 * g(1) + t155 * g(2);
t125 = t151 * g(1) + t154 * g(2);
t124 = t150 * g(1) + t153 * g(2);
t105 = t158 - t123;
t104 = t157 - t122;
t103 = t156 - t121;
t72 = -t327 - t370;
t71 = -t328 - t371;
t70 = -t329 - t372;
t66 = t91 * t198 + t242 * t367;
t65 = t93 * t198 + t240 * t369;
t64 = t92 * t198 + t241 * t368;
t63 = -t187 * t344 + t75;
t62 = -t185 * t347 + t74;
t61 = -t183 * t350 + t73;
t60 = -t327 + (-t262 - t344) * t193;
t59 = -t328 + (-t264 - t347) * t191;
t58 = -t329 + (-t266 - t350) * t189;
t30 = t149 * t240 - t187 * t330;
t29 = t148 * t241 - t185 * t331;
t28 = t147 * t242 - t183 * t332;
t24 = t198 * t262 + t27;
t23 = t198 * t264 + t26;
t22 = t198 * t266 + t25;
t21 = t187 * t65 + t24 * t193;
t20 = t24 * t187 - t65 * t193;
t19 = t185 * t64 + t23 * t191;
t18 = t23 * t185 - t64 * t191;
t17 = t183 * t66 + t22 * t189;
t16 = t22 * t183 - t66 * t189;
t15 = t123 - t227 - 0.2e1 * t363;
t14 = t122 - t228 - 0.2e1 * t364;
t13 = t121 - t229 - 0.2e1 * t365;
t12 = -t226 - t378;
t11 = -t225 - t377;
t10 = -t224 - t376;
t9 = t126 * t193 - t379 * t187;
t8 = t126 * t187 + t379 * t193;
t7 = t125 * t191 - t380 * t185;
t6 = t125 * t185 + t380 * t191;
t5 = t124 * t189 - t381 * t183;
t4 = t124 * t183 + t381 * t189;
t3 = t199 * t188 + t39 * t202 + t382 * qJ(2,1) + (t234 + t378) * pkin(1);
t2 = t199 * t186 + t38 * t201 + t383 * qJ(2,2) + (t235 + t377) * pkin(1);
t1 = t199 * t184 + t37 * t200 + t384 * qJ(2,3) + (t236 + t376) * pkin(1);
t31 = [t39 * t255 + t38 * t257 + t37 * t259, t103 * t259 + t104 * t257 + t105 * t255, t106 * t259 + t107 * t257 + t108 * t255, (t15 * t294 + t39 * t334) * t139 + (t14 * t296 + t38 * t336) * t136 + (t13 * t298 + t37 * t338) * t133, t25 * t259 + t27 * t255 + t26 * t257 - t82 * t273 - t84 * t276 - t83 * t279, (t12 * t334 + t3 * t294) * t139 + (t11 * t336 + t2 * t296) * t136 + (t1 * t298 + t10 * t338) * t133, t153 * t249 + t154 * t247 + t155 * t245 + (-t150 * t280 - t151 * t277 - t152 * t274) * t206, (-t150 * t248 - t151 * t246 - t152 * t244) * t206 + 0.2e1 * t28 * t259 + 0.2e1 * t29 * t257 + 0.2e1 * t30 * t255, t73 * t259 + t74 * t257 + t75 * t255 + (-t150 * t272 - t151 * t271 - t152 * t270) * t206, t70 * t259 + t71 * t257 + t72 * t255 + (t150 * t37 + t151 * t38 + t152 * t39) * t206, (-t93 * t301 - t92 * t303 - t91 * t305) * t206, (t20 * t294 + t63 * t334) * t139 + (t18 * t296 + t62 * t336) * t136 + (t16 * t298 + t61 * t338) * t133 + (-t8 * t301 - t6 * t303 - t4 * t305) * t206, (t21 * t294 + t60 * t334) * t139 + (t19 * t296 + t59 * t336) * t136 + (t17 * t298 + t58 * t338) * t133 + (-t9 * t301 - t7 * t303 - t5 * t305) * t206, t182 - g(1); -t39 * t256 - t38 * t258 - t37 * t260, -t103 * t260 - t104 * t258 - t105 * t256, -t106 * t260 - t107 * t258 - t108 * t256, (-t15 * t300 + t39 * t333) * t139 + (-t14 * t302 + t38 * t335) * t136 + (-t13 * t304 + t37 * t337) * t133, -t25 * t260 - t27 * t256 - t26 * t258 - t85 * t273 - t87 * t276 - t86 * t279, (t12 * t333 - t3 * t300) * t139 + (t11 * t335 - t2 * t302) * t136 + (-t1 * t304 + t10 * t337) * t133, -t150 * t249 - t151 * t247 - t152 * t245 + (-t153 * t280 - t154 * t277 - t155 * t274) * t206, (-t153 * t248 - t154 * t246 - t155 * t244) * t206 - 0.2e1 * t28 * t260 - 0.2e1 * t29 * t258 - 0.2e1 * t30 * t256, -t73 * t260 - t74 * t258 - t75 * t256 + (-t153 * t272 - t154 * t271 - t155 * t270) * t206, -t70 * t260 - t71 * t258 - t72 * t256 + (t153 * t37 + t154 * t38 + t155 * t39) * t206, (-t93 * t295 - t92 * t297 - t91 * t299) * t206, (-t20 * t300 + t63 * t333) * t139 + (-t18 * t302 + t62 * t335) * t136 + (-t16 * t304 + t61 * t337) * t133 + (-t8 * t295 - t6 * t297 - t4 * t299) * t206, (-t21 * t300 + t60 * t333) * t139 + (-t19 * t302 + t59 * t335) * t136 + (-t17 * t304 + t58 * t337) * t133 + (-t295 * t9 - t297 * t7 - t299 * t5) * t206, t181 - g(2); t39 * t307 + t38 * t311 + t37 * t315, t103 * t315 + t104 * t311 + t105 * t307, t106 * t315 + t107 * t311 + t108 * t307, (t111 * t39 + t15 * t194) * t139 + (t110 * t38 + t14 * t192) * t136 + (t109 * t37 + t13 * t190) * t133, -t109 * t348 - t110 * t345 - t111 * t342 + t25 * t315 + t26 * t311 + t27 * t307, (t111 * t12 + t194 * t3) * t139 + (t11 * t110 + t192 * t2) * t136 + (t1 * t190 + t10 * t109) * t133, t190 * t281 + t192 * t278 + t194 * t275, 0.2e1 * t28 * t315 + 0.2e1 * t29 * t311 + 0.2e1 * t30 * t307, t75 * t307 + t74 * t311 + t73 * t315, t72 * t307 + t71 * t311 + t70 * t315, 0, (t111 * t63 + t194 * t20) * t139 + (t110 * t62 + t18 * t192) * t136 + (t109 * t61 + t16 * t190) * t133, (t111 * t60 + t194 * t21) * t139 + (t110 * t59 + t19 * t192) * t136 + (t109 * t58 + t17 * t190) * t133, t180 - g(3);];
tauX_reg  = t31;
