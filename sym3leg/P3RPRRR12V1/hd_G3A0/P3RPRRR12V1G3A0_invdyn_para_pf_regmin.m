% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR12V1G3A0
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR12V1G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:13
% EndTime: 2020-08-06 18:28:19
% DurationCPUTime: 5.38s
% Computational Cost: add. (15559->390), mult. (25515->732), div. (3228->14), fcn. (19026->18), ass. (0->299)
t173 = legFrame(3,2);
t143 = sin(t173);
t146 = cos(t173);
t124 = t146 * g(1) - t143 * g(2);
t180 = sin(qJ(1,3));
t186 = cos(qJ(1,3));
t104 = -g(3) * t180 + t124 * t186;
t179 = sin(qJ(3,3));
t356 = t179 * pkin(3);
t136 = qJ(2,3) + t356;
t128 = 0.1e1 / t136 ^ 2;
t163 = 0.1e1 / t179;
t317 = t128 * t163;
t191 = xDP(3);
t155 = pkin(3) * t191;
t360 = (-pkin(5) - pkin(6));
t162 = pkin(1) - t360;
t192 = xDP(2);
t193 = xDP(1);
t235 = -t143 * t192 + t146 * t193;
t185 = cos(qJ(3,3));
t253 = t180 * pkin(3) * (t185 - 0.1e1) * (t185 + 0.1e1);
t115 = t143 * t193 + t146 * t192;
t285 = t185 * t115;
t293 = t162 * t191;
t169 = t185 ^ 2;
t375 = -t169 + 0.1e1;
t52 = ((qJ(2,3) * t191 + t235 * t162) * t186 + (t235 * qJ(2,3) - t293) * t180 + pkin(3) * t285) * t179 - t235 * t253 + qJ(2,3) * t285 + t375 * t186 * t155;
t97 = -t180 * t191 + t235 * t186;
t378 = 0.2e1 * t97 * t52 * t317 - t104;
t174 = legFrame(2,2);
t144 = sin(t174);
t147 = cos(t174);
t125 = t147 * g(1) - t144 * g(2);
t182 = sin(qJ(1,2));
t188 = cos(qJ(1,2));
t106 = -g(3) * t182 + t125 * t188;
t181 = sin(qJ(3,2));
t355 = t181 * pkin(3);
t137 = qJ(2,2) + t355;
t131 = 0.1e1 / t137 ^ 2;
t165 = 0.1e1 / t181;
t313 = t131 * t165;
t234 = -t144 * t192 + t147 * t193;
t187 = cos(qJ(3,2));
t252 = t182 * pkin(3) * (t187 - 0.1e1) * (t187 + 0.1e1);
t117 = t144 * t193 + t147 * t192;
t283 = t187 * t117;
t170 = t187 ^ 2;
t374 = -t170 + 0.1e1;
t53 = ((qJ(2,2) * t191 + t234 * t162) * t188 + (t234 * qJ(2,2) - t293) * t182 + pkin(3) * t283) * t181 - t234 * t252 + qJ(2,2) * t283 + t374 * t188 * t155;
t98 = -t182 * t191 + t234 * t188;
t377 = 0.2e1 * t98 * t53 * t313 - t106;
t175 = legFrame(1,2);
t145 = sin(t175);
t148 = cos(t175);
t126 = t148 * g(1) - t145 * g(2);
t184 = sin(qJ(1,1));
t190 = cos(qJ(1,1));
t108 = -g(3) * t184 + t126 * t190;
t183 = sin(qJ(3,1));
t354 = t183 * pkin(3);
t138 = qJ(2,1) + t354;
t134 = 0.1e1 / t138 ^ 2;
t167 = 0.1e1 / t183;
t309 = t134 * t167;
t233 = -t145 * t192 + t148 * t193;
t189 = cos(qJ(3,1));
t251 = t184 * pkin(3) * (t189 - 0.1e1) * (t189 + 0.1e1);
t119 = t145 * t193 + t148 * t192;
t281 = t189 * t119;
t171 = t189 ^ 2;
t373 = -t171 + 0.1e1;
t54 = ((qJ(2,1) * t191 + t233 * t162) * t190 + (t233 * qJ(2,1) - t293) * t184 + pkin(3) * t281) * t183 - t233 * t251 + qJ(2,1) * t281 + t373 * t190 * t155;
t99 = -t184 * t191 + t233 * t190;
t376 = 0.2e1 * t99 * t54 * t309 - t108;
t164 = 0.1e1 / t179 ^ 2;
t203 = 0.1e1 / pkin(3) ^ 2;
t329 = t115 ^ 2 * t203;
t265 = t164 * t329;
t166 = 0.1e1 / t181 ^ 2;
t328 = t117 ^ 2 * t203;
t263 = t166 * t328;
t168 = 0.1e1 / t183 ^ 2;
t327 = t119 ^ 2 * t203;
t261 = t168 * t327;
t372 = t185 * t265;
t371 = t187 * t263;
t370 = t189 * t261;
t226 = t126 * t184;
t133 = 0.1e1 / t138;
t135 = t133 * t134;
t176 = xDDP(3);
t177 = xDDP(2);
t178 = xDDP(1);
t342 = t162 * t99;
t346 = t134 * t99;
t39 = (0.2e1 * t134 * t281 - t135 * t54) * t99 * t167 + (-t184 * t176 - (t167 * t54 - t342) * t346 + (-t145 * t177 + t148 * t178) * t190) * t133;
t357 = pkin(1) * t39;
t369 = -t226 - t357;
t194 = pkin(1) + pkin(5);
t154 = g(3) * t190;
t111 = t138 * t190 - t162 * t184;
t198 = qJ(2,1) ^ 2;
t201 = pkin(3) ^ 2;
t229 = -(pkin(6) ^ 2) - t201 + (-(2 * pkin(6)) - pkin(5)) * pkin(5) + ((2 * t360) - pkin(1)) * pkin(1);
t294 = t162 * t167;
t266 = t54 * t294;
t202 = 0.1e1 / pkin(3);
t290 = t167 * t202;
t310 = t133 * t189;
t324 = t119 * t133;
t238 = t184 * qJ(2,1) + t162 * t190;
t282 = t189 * qJ(2,1);
t287 = t183 * t189;
t84 = t238 * t148 * t183 + t145 * t282 + (t373 * t148 * t184 + t145 * t287) * pkin(3);
t87 = (t148 * pkin(3) * t189 - t238 * t145) * t183 + t145 * t251 + t148 * t282;
t230 = -t281 * t294 * t346 + t99 * t135 * t266 - ((t119 * t167 + t310 * t342) * t183 + qJ(2,1) * t119 * t290) * t168 * t324 + ((t266 + (-0.2e1 * qJ(2,1) * t354 + t171 * t201 - t198 + t229) * t99) * t346 - t111 * t176 + (-t87 * t177 - t84 * t178) * t167) * t133;
t223 = t154 + t230;
t96 = t99 ^ 2;
t347 = t134 * t96;
t222 = qJ(2,1) * t347 + t223;
t368 = t194 * t39 + t222 + t226;
t227 = t125 * t182;
t130 = 0.1e1 / t137;
t132 = t130 * t131;
t343 = t162 * t98;
t349 = t131 * t98;
t38 = (0.2e1 * t131 * t283 - t132 * t53) * t98 * t165 + (-t182 * t176 - (t165 * t53 - t343) * t349 + (-t144 * t177 + t147 * t178) * t188) * t130;
t358 = pkin(1) * t38;
t367 = -t227 - t358;
t153 = g(3) * t188;
t110 = t137 * t188 - t162 * t182;
t197 = qJ(2,2) ^ 2;
t295 = t162 * t165;
t267 = t53 * t295;
t291 = t165 * t202;
t314 = t130 * t187;
t325 = t117 * t130;
t237 = qJ(2,2) * t182 + t162 * t188;
t284 = t187 * qJ(2,2);
t288 = t181 * t187;
t83 = t237 * t147 * t181 + t144 * t284 + (t374 * t147 * t182 + t144 * t288) * pkin(3);
t86 = (t147 * pkin(3) * t187 - t237 * t144) * t181 + t144 * t252 + t147 * t284;
t231 = -t283 * t295 * t349 + t98 * t132 * t267 - ((t117 * t165 + t314 * t343) * t181 + qJ(2,2) * t117 * t291) * t166 * t325 + ((t267 + (-0.2e1 * qJ(2,2) * t355 + t170 * t201 - t197 + t229) * t98) * t349 - t110 * t176 + (-t86 * t177 - t83 * t178) * t165) * t130;
t224 = t153 + t231;
t95 = t98 ^ 2;
t350 = t131 * t95;
t221 = qJ(2,2) * t350 + t224;
t366 = t194 * t38 + t221 + t227;
t228 = t124 * t180;
t127 = 0.1e1 / t136;
t129 = t127 * t128;
t344 = t162 * t97;
t352 = t128 * t97;
t37 = (0.2e1 * t128 * t285 - t129 * t52) * t97 * t163 + (-t180 * t176 - (t163 * t52 - t344) * t352 + (-t143 * t177 + t146 * t178) * t186) * t127;
t359 = pkin(1) * t37;
t365 = -t228 - t359;
t152 = g(3) * t186;
t109 = t136 * t186 - t162 * t180;
t196 = qJ(2,3) ^ 2;
t296 = t162 * t163;
t268 = t52 * t296;
t292 = t163 * t202;
t318 = t127 * t185;
t326 = t115 * t127;
t236 = qJ(2,3) * t180 + t162 * t186;
t286 = t185 * qJ(2,3);
t289 = t179 * t185;
t82 = t236 * t146 * t179 + t143 * t286 + (t375 * t146 * t180 + t143 * t289) * pkin(3);
t85 = (t146 * pkin(3) * t185 - t236 * t143) * t179 + t143 * t253 + t146 * t286;
t232 = -t285 * t296 * t352 + t97 * t129 * t268 - ((t115 * t163 + t318 * t344) * t179 + qJ(2,3) * t115 * t292) * t164 * t326 + ((t268 + (-0.2e1 * qJ(2,3) * t356 + t169 * t201 - t196 + t229) * t97) * t352 - t109 * t176 + (-t85 * t177 - t82 * t178) * t163) * t127;
t225 = t152 + t232;
t94 = t97 ^ 2;
t353 = t128 * t94;
t220 = qJ(2,3) * t353 + t225;
t364 = t194 * t37 + t220 + t228;
t363 = 0.2e1 * qJ(2,1);
t362 = 0.2e1 * qJ(2,2);
t361 = 0.2e1 * qJ(2,3);
t140 = 0.2e1 * t169 - 0.1e1;
t141 = 0.2e1 * t170 - 0.1e1;
t142 = 0.2e1 * t171 - 0.1e1;
t351 = t129 * t94;
t348 = t132 * t95;
t345 = t135 * t96;
t341 = t163 * t82;
t340 = t163 * t85;
t339 = t165 * t83;
t338 = t165 * t86;
t337 = t167 * t84;
t336 = t167 * t87;
t335 = t185 * t37;
t334 = t187 * t38;
t333 = t189 * t39;
t91 = -t163 * t372 + (-t143 * t178 - t146 * t177) * t292;
t332 = t91 * t179;
t92 = -t165 * t371 + (-t144 * t178 - t147 * t177) * t291;
t331 = t92 * t181;
t93 = -t167 * t370 + (-t145 * t178 - t148 * t177) * t290;
t330 = t93 * t183;
t319 = t127 * t180;
t315 = t130 * t182;
t311 = t133 * t184;
t308 = t143 * t163;
t307 = t143 * t186;
t306 = t144 * t165;
t305 = t144 * t188;
t304 = t145 * t167;
t303 = t145 * t190;
t302 = t146 * t163;
t301 = t146 * t186;
t300 = t147 * t165;
t299 = t147 * t188;
t298 = t148 * t167;
t297 = t148 * t190;
t250 = t202 * t97 * t326;
t280 = (0.2e1 * t250 + t335) * t318;
t279 = t185 * t353;
t278 = t163 * t351;
t249 = t202 * t98 * t325;
t277 = (0.2e1 * t249 + t334) * t314;
t276 = t187 * t350;
t275 = t165 * t348;
t242 = t99 * t202 * t324;
t274 = (0.2e1 * t242 + t333) * t310;
t273 = t189 * t347;
t272 = t167 * t345;
t271 = t163 * t335;
t270 = t165 * t334;
t269 = t167 * t333;
t259 = t127 * t307;
t258 = t127 * t301;
t257 = t130 * t305;
t256 = t130 * t299;
t255 = t133 * t303;
t254 = t133 * t297;
t248 = t186 * t280;
t247 = t140 * t94 * t317;
t246 = t188 * t277;
t245 = t141 * t95 * t313;
t244 = t190 * t274;
t243 = t142 * t96 * t309;
t73 = -t163 * t329 + t91 * t185;
t74 = -t165 * t328 + t92 * t187;
t75 = -t167 * t327 + t93 * t189;
t241 = t165 * t249;
t240 = t163 * t250;
t239 = t167 * t242;
t25 = t37 * t361 + t378;
t26 = t38 * t362 + t377;
t27 = t39 * t363 + t376;
t195 = pkin(1) * g(3);
t123 = t145 * g(1) + t148 * g(2);
t122 = t144 * g(1) + t147 * g(2);
t121 = t143 * g(1) + t146 * g(2);
t107 = t154 + t226;
t105 = t153 + t227;
t103 = t152 + t228;
t72 = -t330 - t370;
t71 = -t331 - t371;
t70 = -t332 - t372;
t66 = t93 * t194 + t239 * t363;
t65 = t92 * t194 + t241 * t362;
t64 = t91 * t194 + t240 * t361;
t63 = -t183 * t347 + t75;
t62 = -t181 * t350 + t74;
t61 = -t179 * t353 + t73;
t60 = -t330 + (-t261 - t347) * t189;
t59 = -t331 + (-t263 - t350) * t187;
t58 = -t332 + (-t265 - t353) * t185;
t30 = t142 * t239 - t39 * t287;
t29 = t141 * t241 - t38 * t288;
t28 = t140 * t240 - t37 * t289;
t24 = t194 * t261 + t27;
t23 = t194 * t263 + t26;
t22 = t194 * t265 + t25;
t21 = t183 * t66 + t24 * t189;
t20 = t24 * t183 - t189 * t66;
t19 = t181 * t65 + t23 * t187;
t18 = t23 * t181 - t187 * t65;
t17 = t179 * t64 + t22 * t185;
t16 = t22 * t179 - t185 * t64;
t15 = -t223 - t226 - 0.2e1 * t357;
t14 = -t224 - t227 - 0.2e1 * t358;
t13 = -t225 - t228 - 0.2e1 * t359;
t12 = -t222 + t369;
t11 = -t221 + t367;
t10 = -t220 + t365;
t9 = t123 * t183 - t368 * t189;
t8 = t123 * t189 + t368 * t183;
t7 = t122 * t181 - t366 * t187;
t6 = t122 * t187 + t366 * t181;
t5 = t121 * t179 - t364 * t185;
t4 = t121 * t185 + t364 * t179;
t3 = t195 * t190 + t39 * t198 + t376 * qJ(2,1) + (t230 - t369) * pkin(1);
t2 = t195 * t188 + t38 * t197 + t377 * qJ(2,2) + (t231 - t367) * pkin(1);
t1 = t195 * t186 + t37 * t196 + t378 * qJ(2,3) + (t232 - t365) * pkin(1);
t31 = [t39 * t254 + t38 * t256 + t37 * t258, t103 * t258 + t105 * t256 + t107 * t254, t104 * t258 + t106 * t256 + t108 * t254, (t15 * t297 + t39 * t337) * t133 + (t14 * t299 + t38 * t339) * t130 + (t13 * t301 + t37 * t341) * t127, t25 * t258 + t27 * t254 + t26 * t256 - t84 * t272 - t83 * t275 - t82 * t278, (t12 * t337 + t3 * t297) * t133 + (t11 * t339 + t2 * t299) * t130 + (t1 * t301 + t10 * t341) * t127, t146 * t248 + t147 * t246 + t148 * t244 + (-t143 * t279 - t144 * t276 - t145 * t273) * t202, (-t143 * t247 - t144 * t245 - t145 * t243) * t202 + 0.2e1 * t28 * t258 + 0.2e1 * t29 * t256 + 0.2e1 * t30 * t254, t73 * t258 + t74 * t256 + t75 * t254 + (-t143 * t271 - t144 * t270 - t145 * t269) * t202, t70 * t258 + t71 * t256 + t72 * t254 + (t143 * t37 + t144 * t38 + t145 * t39) * t202, (-t93 * t304 - t92 * t306 - t91 * t308) * t202, (t20 * t297 + t63 * t337) * t133 + (t18 * t299 + t62 * t339) * t130 + (t16 * t301 + t61 * t341) * t127 + (-t9 * t304 - t7 * t306 - t5 * t308) * t202, (t21 * t297 + t60 * t337) * t133 + (t19 * t299 + t59 * t339) * t130 + (t17 * t301 + t58 * t341) * t127 + (-t304 * t8 - t306 * t6 - t308 * t4) * t202, t178 - g(1); -t39 * t255 - t38 * t257 - t37 * t259, -t103 * t259 - t105 * t257 - t107 * t255, -t104 * t259 - t106 * t257 - t108 * t255, (-t15 * t303 + t39 * t336) * t133 + (-t14 * t305 + t38 * t338) * t130 + (-t13 * t307 + t37 * t340) * t127, -t25 * t259 - t27 * t255 - t26 * t257 - t87 * t272 - t86 * t275 - t85 * t278, (t12 * t336 - t3 * t303) * t133 + (t11 * t338 - t2 * t305) * t130 + (-t1 * t307 + t10 * t340) * t127, -t143 * t248 - t144 * t246 - t145 * t244 + (-t146 * t279 - t147 * t276 - t148 * t273) * t202, (-t146 * t247 - t147 * t245 - t148 * t243) * t202 - 0.2e1 * t28 * t259 - 0.2e1 * t29 * t257 - 0.2e1 * t30 * t255, -t73 * t259 - t74 * t257 - t75 * t255 + (-t146 * t271 - t147 * t270 - t148 * t269) * t202, -t70 * t259 - t71 * t257 - t72 * t255 + (t146 * t37 + t147 * t38 + t148 * t39) * t202, (-t93 * t298 - t92 * t300 - t91 * t302) * t202, (-t20 * t303 + t63 * t336) * t133 + (-t18 * t305 + t62 * t338) * t130 + (-t16 * t307 + t61 * t340) * t127 + (-t9 * t298 - t7 * t300 - t5 * t302) * t202, (-t21 * t303 + t60 * t336) * t133 + (-t19 * t305 + t59 * t338) * t130 + (-t17 * t307 + t58 * t340) * t127 + (-t298 * t8 - t300 * t6 - t302 * t4) * t202, t177 - g(2); -t39 * t311 - t38 * t315 - t37 * t319, -t103 * t319 - t105 * t315 - t107 * t311, -t104 * t319 - t106 * t315 - t108 * t311, (t111 * t39 - t15 * t184) * t133 + (t110 * t38 - t14 * t182) * t130 + (t109 * t37 - t13 * t180) * t127, -t109 * t351 - t110 * t348 - t111 * t345 - t25 * t319 - t26 * t315 - t27 * t311, (t111 * t12 - t184 * t3) * t133 + (t11 * t110 - t182 * t2) * t130 + (-t1 * t180 + t10 * t109) * t127, -t180 * t280 - t182 * t277 - t184 * t274, -0.2e1 * t28 * t319 - 0.2e1 * t29 * t315 - 0.2e1 * t30 * t311, -t75 * t311 - t74 * t315 - t73 * t319, -t72 * t311 - t71 * t315 - t70 * t319, 0, (t111 * t63 - t184 * t20) * t133 + (t110 * t62 - t18 * t182) * t130 + (t109 * t61 - t16 * t180) * t127, (t111 * t60 - t184 * t21) * t133 + (t110 * t59 - t182 * t19) * t130 + (t109 * t58 - t17 * t180) * t127, t176 - g(3);];
tauX_reg  = t31;
