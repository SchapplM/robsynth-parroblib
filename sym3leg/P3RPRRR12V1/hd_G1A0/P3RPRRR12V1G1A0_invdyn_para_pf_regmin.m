% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR12V1G1A0
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
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR12V1G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:30
% EndTime: 2020-08-06 18:21:34
% DurationCPUTime: 3.59s
% Computational Cost: add. (12023->349), mult. (18258->627), div. (2190->14), fcn. (14151->18), ass. (0->263)
t187 = sin(qJ(3,3));
t172 = 0.1e1 / t187 ^ 2;
t199 = xDP(3);
t272 = t199 ^ 2 / pkin(3) ^ 2;
t332 = t172 * t272;
t189 = sin(qJ(3,2));
t174 = 0.1e1 / t189 ^ 2;
t331 = t174 * t272;
t191 = sin(qJ(3,1));
t176 = 0.1e1 / t191 ^ 2;
t330 = t176 * t272;
t193 = cos(qJ(3,3));
t329 = t193 * t332;
t195 = cos(qJ(3,2));
t328 = t195 * t331;
t197 = cos(qJ(3,1));
t327 = t197 * t330;
t183 = legFrame(1,3);
t156 = sin(t183);
t159 = cos(t183);
t117 = t156 * g(1) - t159 * g(2);
t120 = t159 * g(1) + t156 * g(2);
t192 = sin(qJ(1,1));
t198 = cos(qJ(1,1));
t87 = -t117 * t192 + t120 * t198;
t182 = legFrame(2,3);
t155 = sin(t182);
t158 = cos(t182);
t116 = t155 * g(1) - t158 * g(2);
t119 = t158 * g(1) + t155 * g(2);
t190 = sin(qJ(1,2));
t196 = cos(qJ(1,2));
t86 = -t116 * t190 + t119 * t196;
t181 = legFrame(3,3);
t154 = sin(t181);
t157 = cos(t181);
t115 = t154 * g(1) - t157 * g(2);
t118 = t157 * g(1) + t154 * g(2);
t188 = sin(qJ(1,3));
t194 = cos(qJ(1,3));
t85 = -t115 * t188 + t118 * t194;
t150 = t187 * pkin(3) + qJ(2,3);
t139 = 0.1e1 / t150;
t151 = t189 * pkin(3) + qJ(2,2);
t142 = 0.1e1 / t151;
t152 = -t191 * pkin(3) - qJ(2,1);
t145 = 0.1e1 / t152;
t171 = 0.1e1 / t187;
t173 = 0.1e1 / t189;
t175 = 0.1e1 / t191;
t140 = 0.1e1 / t150 ^ 2;
t143 = 0.1e1 / t151 ^ 2;
t146 = 0.1e1 / t152 ^ 2;
t326 = 0.2e1 * qJ(2,1);
t325 = 0.2e1 * qJ(2,2);
t324 = 0.2e1 * qJ(2,3);
t202 = pkin(1) + pkin(5);
t278 = t157 * t194;
t106 = -t154 * t188 + t278;
t281 = t154 * t194;
t109 = t157 * t188 + t281;
t141 = t139 * t140;
t185 = xDDP(2);
t186 = xDDP(1);
t267 = t199 * t193;
t170 = pkin(6) + t202;
t201 = xDP(1);
t148 = t170 * t201;
t200 = xDP(2);
t121 = -qJ(2,3) * t200 + t148;
t149 = t200 * t170;
t124 = qJ(2,3) * t201 + t149;
t236 = t188 * t200 + t194 * t201;
t237 = -t188 * t201 + t194 * t200;
t46 = ((t121 * t194 + t188 * t124) * t157 + (-t188 * t121 + t124 * t194) * t154) * t187 + qJ(2,3) * t267 + (t187 * t267 + (-t236 * t154 + t237 * t157) * (t193 - 0.1e1) * (t193 + 0.1e1)) * pkin(3);
t305 = t171 * t46;
t79 = t237 * t154 + t236 * t157;
t308 = t170 * t79;
t31 = (0.2e1 * t140 * t267 - t141 * t46) * t79 * t171 + (t106 * t186 + t109 * t185 - (t305 - t308) * t140 * t79) * t139;
t323 = pkin(1) * t31;
t277 = t158 * t196;
t107 = -t155 * t190 + t277;
t280 = t155 * t196;
t110 = t158 * t190 + t280;
t144 = t142 * t143;
t266 = t199 * t195;
t122 = -qJ(2,2) * t200 + t148;
t125 = qJ(2,2) * t201 + t149;
t234 = t190 * t200 + t196 * t201;
t235 = -t190 * t201 + t196 * t200;
t47 = ((t122 * t196 + t190 * t125) * t158 + (-t190 * t122 + t125 * t196) * t155) * t189 + qJ(2,2) * t266 + (t189 * t266 + (-t234 * t155 + t235 * t158) * (t195 - 0.1e1) * (t195 + 0.1e1)) * pkin(3);
t304 = t173 * t47;
t80 = t235 * t155 + t234 * t158;
t307 = t170 * t80;
t32 = (0.2e1 * t143 * t266 - t144 * t47) * t80 * t173 + (t107 * t186 + t110 * t185 - (t304 - t307) * t143 * t80) * t142;
t322 = pkin(1) * t32;
t276 = t159 * t198;
t108 = -t156 * t192 + t276;
t279 = t156 * t198;
t111 = t159 * t192 + t279;
t147 = t145 * t146;
t265 = t199 * t197;
t123 = -qJ(2,1) * t200 + t148;
t126 = qJ(2,1) * t201 + t149;
t232 = t192 * t200 + t198 * t201;
t233 = -t192 * t201 + t198 * t200;
t48 = ((t123 * t198 + t192 * t126) * t159 + (-t192 * t123 + t126 * t198) * t156) * t191 + qJ(2,1) * t265 + (t191 * t265 + (-t232 * t156 + t233 * t159) * (t197 - 0.1e1) * (t197 + 0.1e1)) * pkin(3);
t303 = t175 * t48;
t81 = t233 * t156 + t232 * t159;
t306 = t170 * t81;
t33 = (0.2e1 * t146 * t265 + t147 * t48) * t81 * t175 - (t108 * t186 + t111 * t185 - (t303 - t306) * t146 * t81) * t145;
t321 = pkin(1) * t33;
t177 = t193 ^ 2;
t320 = 0.2e1 * t177 - 0.1e1;
t178 = t195 ^ 2;
t319 = 0.2e1 * t178 - 0.1e1;
t179 = t197 ^ 2;
t318 = 0.2e1 * t179 - 0.1e1;
t286 = t140 * t171;
t261 = t79 * t286;
t43 = 0.2e1 * t46 * t261;
t317 = t31 * t324 + t43;
t284 = t143 * t173;
t259 = t80 * t284;
t44 = 0.2e1 * t47 * t259;
t316 = t32 * t325 + t44;
t282 = t146 * t175;
t257 = t81 * t282;
t45 = 0.2e1 * t48 * t257;
t315 = t33 * t326 + t45;
t76 = t79 ^ 2;
t314 = t140 * t76;
t313 = t141 * t76;
t77 = t80 ^ 2;
t312 = t143 * t77;
t311 = t144 * t77;
t78 = t81 ^ 2;
t310 = t146 * t78;
t309 = t147 * t78;
t184 = xDDP(3);
t207 = 0.1e1 / pkin(3);
t271 = t184 * t207;
t97 = (-t271 - t329) * t171;
t302 = t97 * t187;
t98 = (-t271 - t328) * t173;
t301 = t98 * t189;
t99 = (-t271 - t327) * t175;
t300 = t99 * t191;
t299 = t106 * t139;
t298 = t107 * t142;
t297 = t108 * t145;
t296 = t109 * t139;
t295 = t110 * t142;
t294 = t111 * t145;
t112 = t118 * t188;
t113 = t119 * t190;
t114 = t120 * t192;
t287 = t139 * t193;
t285 = t142 * t195;
t283 = t145 * t197;
t275 = t171 * t193;
t274 = t173 * t195;
t273 = t175 * t197;
t270 = t187 * t193;
t269 = t189 * t195;
t268 = t191 * t197;
t264 = t199 * t207;
t206 = pkin(3) ^ 2;
t263 = -t170 ^ 2 - t206;
t243 = t79 * t139 * t264;
t262 = (t193 * t31 + 0.2e1 * t243) * t287;
t242 = t80 * t142 * t264;
t260 = (t195 * t32 + 0.2e1 * t242) * t285;
t241 = t81 * t145 * t264;
t258 = (t197 * t33 - 0.2e1 * t241) * t283;
t256 = t170 * t305;
t255 = t170 * t304;
t254 = t170 * t303;
t253 = qJ(2,3) * t314;
t252 = t193 * t314;
t251 = qJ(2,2) * t312;
t250 = t195 * t312;
t249 = qJ(2,1) * t310;
t248 = t197 * t310;
t244 = t202 * t272;
t91 = -t171 * t272 + t97 * t193;
t92 = -t173 * t272 + t98 * t195;
t93 = -t175 * t272 + t99 * t197;
t240 = t175 * t241;
t239 = t171 * t243;
t238 = t173 * t242;
t203 = qJ(2,3) ^ 2;
t228 = pkin(3) * t170 * t264;
t100 = t188 * t150 + t170 * t194;
t103 = -t150 * t194 + t170 * t188;
t70 = t100 * t157 - t154 * t103;
t73 = t154 * t100 + t103 * t157;
t231 = t171 * t228 * t270 * t261 - t79 * t141 * t256 + t184 * t275 + (((t171 * t199 + t287 * t308) * t187 + qJ(2,3) * t171 * t264) * t172 * t199 + t70 * t186 + t73 * t185 - (t187 * t256 + ((t177 * t206 - t203 + t263) * t187 + (t177 - 0.1e1) * t324 * pkin(3)) * t79) * t261) * t139;
t204 = qJ(2,2) ^ 2;
t101 = t190 * t151 + t170 * t196;
t104 = -t151 * t196 + t170 * t190;
t71 = t101 * t158 - t155 * t104;
t74 = t155 * t101 + t104 * t158;
t230 = t173 * t228 * t269 * t259 - t80 * t144 * t255 + t184 * t274 + (((t173 * t199 + t285 * t307) * t189 + qJ(2,2) * t173 * t264) * t174 * t199 + t71 * t186 + t74 * t185 - (t189 * t255 + ((t178 * t206 - t204 + t263) * t189 + (t178 - 0.1e1) * t325 * pkin(3)) * t80) * t259) * t142;
t205 = qJ(2,1) ^ 2;
t102 = -t192 * t152 + t170 * t198;
t105 = t152 * t198 + t170 * t192;
t72 = t102 * t159 - t156 * t105;
t75 = t156 * t102 + t105 * t159;
t229 = t175 * t228 * t268 * t257 + t81 * t147 * t254 + t184 * t273 + (-((t175 * t199 - t283 * t306) * t191 + qJ(2,1) * t175 * t264) * t176 * t199 - t72 * t186 - t75 * t185 - (-t191 * t254 - ((t179 * t206 - t205 + t263) * t191 + (t179 - 0.1e1) * t326 * pkin(3)) * t81) * t257) * t145;
t227 = t33 * t273 + t32 * t274 + t31 * t275;
t226 = -g(1) * t281 + g(2) * t278 - t112 + t231;
t225 = -g(1) * t280 + g(2) * t277 - t113 + t230;
t224 = -g(1) * t279 + g(2) * t276 - t114 + t229;
t131 = g(1) * t192 - g(2) * t198;
t132 = g(1) * t198 + g(2) * t192;
t223 = t131 * t159 + t132 * t156 + t202 * t33 - t229 + t249;
t129 = -t190 * g(1) + t196 * g(2);
t130 = t196 * g(1) + t190 * g(2);
t222 = -t129 * t158 + t130 * t155 + t202 * t32 - t230 + t251;
t127 = -t188 * g(1) + t194 * g(2);
t128 = t194 * g(1) + t188 * g(2);
t221 = -t127 * t157 + t128 * t154 + t202 * t31 - t231 + t253;
t90 = -t300 - t327;
t89 = -t301 - t328;
t88 = -t302 - t329;
t84 = t117 * t198 + t114;
t83 = t116 * t196 + t113;
t82 = t115 * t194 + t112;
t60 = -0.2e1 * qJ(2,1) * t240 + t99 * t202;
t59 = t98 * t202 + t238 * t325;
t58 = t97 * t202 + t239 * t324;
t57 = -t191 * t310 + t93;
t56 = -t189 * t312 + t92;
t55 = -t187 * t314 + t91;
t54 = -t300 + (-t310 - t330) * t197;
t53 = -t301 + (-t312 - t331) * t195;
t52 = -t302 + (-t314 - t332) * t193;
t24 = -t318 * t240 - t33 * t268;
t23 = t319 * t238 - t32 * t269;
t22 = t320 * t239 - t31 * t270;
t21 = t315 - t87;
t20 = t316 - t86;
t19 = t317 - t85;
t18 = t131 * t156 - t132 * t159 + t176 * t244 + t315;
t17 = -t129 * t155 - t130 * t158 + t174 * t244 + t316;
t16 = -t127 * t154 - t128 * t157 + t172 * t244 + t317;
t15 = t18 * t197 + t191 * t60;
t14 = t18 * t191 - t197 * t60;
t13 = t17 * t195 + t189 * t59;
t12 = t17 * t189 - t195 * t59;
t11 = t16 * t193 + t187 * t58;
t10 = t16 * t187 - t193 * t58;
t9 = t224 - 0.2e1 * t321;
t8 = t225 - 0.2e1 * t322;
t7 = t226 - 0.2e1 * t323;
t6 = t224 - t249 - t321;
t5 = t225 - t251 - t322;
t4 = t226 - t253 - t323;
t3 = t33 * t205 + qJ(2,1) * t45 + (-t229 + t321) * pkin(1) + t120 * (t192 * pkin(1) - t198 * qJ(2,1)) + t117 * (t198 * pkin(1) + t192 * qJ(2,1));
t2 = t32 * t204 + qJ(2,2) * t44 + (-t230 + t322) * pkin(1) + t119 * (t190 * pkin(1) - t196 * qJ(2,2)) + t116 * (t196 * pkin(1) + t190 * qJ(2,2));
t1 = t31 * t203 + qJ(2,3) * t43 + (-t231 + t323) * pkin(1) + t118 * (t188 * pkin(1) - t194 * qJ(2,3)) + t115 * (t194 * pkin(1) + t188 * qJ(2,3));
t25 = [-t33 * t297 + t32 * t298 + t31 * t299, -t84 * t297 + t83 * t298 + t82 * t299, -t87 * t297 + t86 * t298 + t85 * t299, -(t108 * t9 + t33 * t72) * t145 + (t107 * t8 + t32 * t71) * t142 + (t106 * t7 + t31 * t70) * t139, t19 * t299 + t20 * t298 - t21 * t297 + t72 * t309 - t71 * t311 - t70 * t313, -(t108 * t3 + t6 * t72) * t145 + (t107 * t2 + t5 * t71) * t142 + (t1 * t106 + t4 * t70) * t139, t106 * t262 + t107 * t260 - t108 * t258, 0.2e1 * t22 * t299 + 0.2e1 * t23 * t298 - 0.2e1 * t24 * t297, -t93 * t297 + t92 * t298 + t91 * t299, -t90 * t297 + t89 * t298 + t88 * t299, 0, -(t108 * t14 + t57 * t72) * t145 + (t107 * t12 + t56 * t71) * t142 + (t10 * t106 + t55 * t70) * t139, -(t108 * t15 + t54 * t72) * t145 + (t107 * t13 + t53 * t71) * t142 + (t106 * t11 + t52 * t70) * t139, t186 - g(1); -t33 * t294 + t32 * t295 + t31 * t296, -t84 * t294 + t83 * t295 + t82 * t296, -t87 * t294 + t86 * t295 + t85 * t296, -(t111 * t9 + t33 * t75) * t145 + (t110 * t8 + t32 * t74) * t142 + (t109 * t7 + t31 * t73) * t139, t19 * t296 + t20 * t295 - t21 * t294 + t75 * t309 - t74 * t311 - t73 * t313, -(t111 * t3 + t6 * t75) * t145 + (t110 * t2 + t5 * t74) * t142 + (t1 * t109 + t4 * t73) * t139, t109 * t262 + t110 * t260 - t111 * t258, 0.2e1 * t22 * t296 + 0.2e1 * t23 * t295 - 0.2e1 * t24 * t294, -t93 * t294 + t92 * t295 + t91 * t296, -t90 * t294 + t89 * t295 + t88 * t296, 0, -(t111 * t14 + t57 * t75) * t145 + (t110 * t12 + t56 * t74) * t142 + (t10 * t109 + t55 * t73) * t139, -(t111 * t15 + t54 * t75) * t145 + (t110 * t13 + t53 * t74) * t142 + (t109 * t11 + t52 * t73) * t139, t185 - g(2); 0, 0, 0, t227, -t171 * t252 - t173 * t250 - t175 * t248, t6 * t273 + t5 * t274 + t4 * t275, (-t248 - t250 - t252) * t207, (-t78 * t318 * t282 - t77 * t319 * t284 - t76 * t320 * t286) * t207, -t227 * t207, (t31 + t32 + t33) * t207, (-t171 * t97 - t173 * t98 - t175 * t99) * t207, t55 * t275 + t56 * t274 + t57 * t273 + (-t175 * (g(3) * t191 - t197 * t223) - t173 * (g(3) * t189 - t195 * t222) - t171 * (g(3) * t187 - t193 * t221)) * t207, t52 * t275 + t53 * t274 + t54 * t273 + (-t175 * (g(3) * t197 + t191 * t223) - t173 * (g(3) * t195 + t189 * t222) - t171 * (g(3) * t193 + t187 * t221)) * t207, t184 - g(3);];
tauX_reg  = t25;
