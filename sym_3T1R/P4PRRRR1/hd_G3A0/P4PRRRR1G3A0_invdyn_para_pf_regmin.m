% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR1G3A0
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
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P4PRRRR1G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:04:55
% EndTime: 2020-03-02 19:05:03
% DurationCPUTime: 8.71s
% Computational Cost: add. (7660->406), mult. (17432->780), div. (4732->22), fcn. (13960->26), ass. (0->300)
t148 = sin(qJ(2,4));
t119 = 0.1e1 / t148 ^ 2;
t151 = legFrame(4,2);
t107 = sin(t151);
t111 = cos(t151);
t147 = sin(qJ(3,4));
t149 = cos(qJ(3,4));
t171 = xDP(4);
t173 = xDP(2);
t174 = xDP(1);
t150 = cos(qJ(2,4));
t172 = xDP(3);
t267 = t150 * t172;
t175 = xP(4);
t115 = sin(t175);
t116 = cos(t175);
t176 = koppelP(4,2);
t180 = koppelP(4,1);
t93 = t115 * t180 + t116 * t176;
t97 = -t115 * t176 + t116 * t180;
t57 = t107 * t97 + t93 * t111;
t50 = (t107 * t173 - t174 * t111 + t57 * t171) * t149 - t147 * t267;
t45 = t50 ^ 2;
t314 = t119 * t45;
t160 = sin(qJ(2,3));
t126 = 0.1e1 / t160 ^ 2;
t152 = legFrame(3,2);
t108 = sin(t152);
t112 = cos(t152);
t159 = sin(qJ(3,3));
t165 = cos(qJ(3,3));
t166 = cos(qJ(2,3));
t265 = t166 * t172;
t177 = koppelP(3,2);
t181 = koppelP(3,1);
t94 = t115 * t181 + t116 * t177;
t98 = -t115 * t177 + t116 * t181;
t58 = t108 * t98 + t94 * t112;
t54 = (t108 * t173 - t174 * t112 + t58 * t171) * t165 - t159 * t265;
t47 = t54 ^ 2;
t310 = t126 * t47;
t162 = sin(qJ(2,2));
t129 = 0.1e1 / t162 ^ 2;
t153 = legFrame(2,2);
t109 = sin(t153);
t113 = cos(t153);
t161 = sin(qJ(3,2));
t167 = cos(qJ(3,2));
t168 = cos(qJ(2,2));
t264 = t168 * t172;
t178 = koppelP(2,2);
t182 = koppelP(2,1);
t95 = t115 * t182 + t116 * t178;
t99 = -t115 * t178 + t116 * t182;
t59 = t109 * t99 + t95 * t113;
t55 = (t109 * t173 - t174 * t113 + t59 * t171) * t167 - t161 * t264;
t48 = t55 ^ 2;
t306 = t129 * t48;
t164 = sin(qJ(2,1));
t132 = 0.1e1 / t164 ^ 2;
t154 = legFrame(1,2);
t110 = sin(t154);
t114 = cos(t154);
t163 = sin(qJ(3,1));
t169 = cos(qJ(3,1));
t170 = cos(qJ(2,1));
t263 = t170 * t172;
t179 = koppelP(1,2);
t183 = koppelP(1,1);
t100 = -t115 * t179 + t116 * t183;
t96 = t115 * t183 + t116 * t179;
t60 = t110 * t100 + t96 * t114;
t56 = (t110 * t173 - t174 * t114 + t60 * t171) * t169 - t163 * t263;
t49 = t56 ^ 2;
t302 = t132 * t49;
t184 = 1 / pkin(2);
t118 = 0.1e1 / t148;
t120 = 0.1e1 / t149;
t125 = 0.1e1 / t160;
t128 = 0.1e1 / t162;
t131 = 0.1e1 / t164;
t133 = 0.1e1 / t165;
t137 = 0.1e1 / t167;
t141 = 0.1e1 / t169;
t185 = 0.1e1 / pkin(2) ^ 2;
t191 = t149 ^ 2;
t121 = 0.1e1 / t191;
t200 = t165 ^ 2;
t134 = 0.1e1 / t200;
t203 = t167 ^ 2;
t138 = 0.1e1 / t203;
t206 = t169 ^ 2;
t142 = 0.1e1 / t206;
t318 = 2 * t184;
t101 = t111 * g(1) - t107 * g(2);
t122 = t120 * t121;
t156 = xDDP(3);
t146 = t172 ^ 2;
t269 = t146 * t184;
t287 = t120 * t147;
t145 = t171 ^ 2;
t155 = xDDP(4);
t157 = xDDP(2);
t61 = -t145 * t93 + t155 * t97 + t157;
t158 = xDDP(1);
t65 = -t145 * t97 - t155 * t93 + t158;
t83 = -t107 * t150 + t111 * t148;
t84 = t148 * t107 + t111 * t150;
t37 = -t107 * g(1) - t111 * g(2) + (t156 * t287 + t83 * t61 + t84 * t65 + (t184 * t314 + t269) * t122) * t118;
t29 = t101 * t148 + t37 * t150;
t317 = t118 * t29;
t316 = t118 * t83;
t315 = t118 * t84;
t102 = t112 * g(1) - t108 * g(2);
t135 = t133 * t134;
t278 = t133 * t159;
t62 = -t145 * t94 + t155 * t98 + t157;
t66 = -t145 * t98 - t155 * t94 + t158;
t85 = -t108 * t166 + t112 * t160;
t86 = t160 * t108 + t112 * t166;
t38 = -t108 * g(1) - t112 * g(2) + (t156 * t278 + t85 * t62 + t86 * t66 + (t184 * t310 + t269) * t135) * t125;
t31 = t102 * t160 + t38 * t166;
t313 = t125 * t31;
t312 = t125 * t85;
t311 = t125 * t86;
t103 = t113 * g(1) - t109 * g(2);
t139 = t137 * t138;
t275 = t137 * t161;
t63 = -t145 * t95 + t155 * t99 + t157;
t67 = -t145 * t99 - t155 * t95 + t158;
t87 = -t109 * t168 + t113 * t162;
t88 = t162 * t109 + t113 * t168;
t39 = -t109 * g(1) - t113 * g(2) + (t156 * t275 + t87 * t63 + t88 * t67 + (t184 * t306 + t269) * t139) * t128;
t32 = t103 * t162 + t39 * t168;
t309 = t128 * t32;
t308 = t128 * t87;
t307 = t128 * t88;
t104 = t114 * g(1) - t110 * g(2);
t143 = t141 * t142;
t272 = t141 * t163;
t64 = t100 * t155 - t145 * t96 + t157;
t68 = -t145 * t100 - t155 * t96 + t158;
t89 = -t110 * t170 + t114 * t164;
t90 = t164 * t110 + t114 * t170;
t40 = -t110 * g(1) - t114 * g(2) + (t156 * t272 + t89 * t64 + t90 * t68 + (t184 * t302 + t269) * t143) * t131;
t33 = t104 * t164 + t40 * t170;
t305 = t131 * t33;
t304 = t131 * t89;
t303 = t131 * t90;
t286 = t121 * t150;
t241 = t147 * t286;
t123 = 0.1e1 / t191 ^ 2;
t285 = t123 * t185;
t242 = t119 * t285;
t262 = t172 * t184;
t288 = t118 * t150;
t25 = -(-t147 * t148 * t172 + t50 * t288) * t50 * t242 + (-t156 * t241 + (t107 * t61 - t111 * t65 - (-t147 * t50 + t267) * t122 * t262) * t120) * t118 * t184;
t301 = t149 * t25;
t277 = t134 * t166;
t231 = t159 * t277;
t136 = 0.1e1 / t200 ^ 2;
t276 = t136 * t185;
t238 = t126 * t276;
t283 = t125 * t166;
t26 = -(-t159 * t160 * t172 + t54 * t283) * t54 * t238 + (-t156 * t231 + (t108 * t62 - t112 * t66 - (-t159 * t54 + t265) * t135 * t262) * t133) * t125 * t184;
t300 = t165 * t26;
t274 = t138 * t168;
t230 = t161 * t274;
t140 = 0.1e1 / t203 ^ 2;
t273 = t140 * t185;
t235 = t129 * t273;
t281 = t128 * t168;
t27 = -(-t161 * t162 * t172 + t55 * t281) * t55 * t235 + (-t156 * t230 + (t109 * t63 - t113 * t67 - (-t161 * t55 + t264) * t139 * t262) * t137) * t128 * t184;
t299 = t167 * t27;
t271 = t142 * t170;
t229 = t163 * t271;
t144 = 0.1e1 / t206 ^ 2;
t270 = t144 * t185;
t232 = t132 * t270;
t279 = t131 * t170;
t28 = -(-t163 * t164 * t172 + t56 * t279) * t56 * t232 + (-t156 * t229 + (t110 * t64 - t114 * t68 - (-t163 * t56 + t263) * t143 * t262) * t141) * t131 * t184;
t298 = t169 * t28;
t268 = t146 * t185;
t228 = t147 * t268;
t266 = t156 * t184;
t79 = t120 * t266 + t122 * t228;
t297 = t79 * t147;
t296 = t79 * t149;
t227 = t159 * t268;
t80 = t133 * t266 + t135 * t227;
t295 = t80 * t159;
t294 = t80 * t165;
t226 = t161 * t268;
t81 = t137 * t266 + t139 * t226;
t293 = t81 * t161;
t292 = t81 * t167;
t225 = t163 * t268;
t82 = t141 * t266 + t143 * t225;
t291 = t82 * t163;
t290 = t82 * t169;
t289 = t118 * t120;
t284 = t125 * t133;
t282 = t128 * t137;
t280 = t131 * t141;
t261 = t172 * t185;
t260 = t57 * t289;
t259 = t123 * t314;
t258 = t58 * t284;
t257 = t136 * t310;
t256 = t59 * t282;
t255 = t140 * t306;
t254 = t60 * t280;
t253 = t144 * t302;
t252 = t107 * t289;
t251 = t108 * t284;
t250 = t109 * t282;
t249 = t110 * t280;
t248 = t111 * t289;
t247 = t112 * t284;
t246 = t113 * t282;
t245 = t114 * t280;
t244 = t118 * t287;
t243 = t118 * t286;
t240 = t125 * t278;
t239 = t125 * t277;
t237 = t128 * t275;
t236 = t128 * t274;
t234 = t131 * t272;
t233 = t131 * t271;
t224 = t29 * t244;
t223 = t118 * t50 * t261;
t222 = t31 * t240;
t221 = t125 * t54 * t261;
t220 = t32 * t237;
t219 = t128 * t55 * t261;
t218 = t33 * t234;
t217 = t131 * t56 * t261;
t216 = t118 * t241;
t215 = t125 * t231;
t214 = t128 * t230;
t213 = t131 * t229;
t212 = 0.2e1 * t121 * t223;
t211 = 0.2e1 * t134 * t221;
t210 = 0.2e1 * t138 * t219;
t209 = 0.2e1 * t142 * t217;
t130 = t163 ^ 2;
t127 = t161 ^ 2;
t124 = t159 ^ 2;
t117 = t147 ^ 2;
t106 = t158 - g(1);
t105 = t157 - g(2);
t92 = -t115 * t155 - t116 * t145;
t91 = -t115 * t145 + t116 * t155;
t78 = t115 * t105 + t116 * t106;
t77 = t116 * t105 - t115 * t106;
t76 = t141 * t268 + t291;
t75 = t137 * t268 + t293;
t74 = t133 * t268 + t295;
t73 = t120 * t268 + t297;
t72 = -t142 * t225 + t290;
t71 = -t138 * t226 + t292;
t70 = -t134 * t227 + t294;
t69 = -t121 * t228 + t296;
t53 = (t100 * t89 - t90 * t96) * t131;
t52 = (t87 * t99 - t88 * t95) * t128;
t51 = (t85 * t98 - t86 * t94) * t125;
t46 = (t83 * t97 - t84 * t93) * t118;
t44 = (t142 * t146 + t253) * t185;
t43 = (t138 * t146 + t255) * t185;
t42 = (t134 * t146 + t257) * t185;
t41 = (t121 * t146 + t259) * t185;
t36 = t104 * t170 - t40 * t164;
t35 = t103 * t168 - t39 * t162;
t34 = t102 * t166 - t38 * t160;
t30 = t101 * t150 - t37 * t148;
t24 = -t49 * t131 * t270 + t170 * t28;
t23 = -t48 * t128 * t273 + t168 * t27;
t22 = -t47 * t125 * t276 + t166 * t26;
t21 = -t49 * t170 * t232 - t28 * t164;
t20 = -t48 * t168 * t235 - t27 * t162;
t19 = -t47 * t166 * t238 - t26 * t160;
t18 = t130 * t28 + t163 * t209;
t17 = t127 * t27 + t161 * t210;
t16 = t124 * t26 + t159 * t211;
t15 = -t45 * t118 * t285 + t150 * t25;
t14 = -t45 * t150 * t242 - t25 * t148;
t13 = t117 * t25 + t147 * t212;
t12 = t163 * t298 + (0.2e1 * t141 - t143) * t217;
t11 = t161 * t299 + (0.2e1 * t137 - t139) * t219;
t10 = t159 * t300 + (0.2e1 * t133 - t135) * t221;
t9 = t147 * t301 + (0.2e1 * t120 - t122) * t223;
t8 = (t44 * t163 - t290) * t164 - t170 * (t163 * t28 + t209);
t7 = (t43 * t161 - t292) * t162 - t168 * (t161 * t27 + t210);
t6 = (t42 * t159 - t294) * t160 - t166 * (t159 * t26 + t211);
t5 = (-t44 * t169 - t291) * t164 + (-0.2e1 * t163 * t143 * t217 + t298) * t170;
t4 = (-t43 * t167 - t293) * t162 + (-0.2e1 * t161 * t139 * t219 + t299) * t168;
t3 = (-t42 * t165 - t295) * t160 + (-0.2e1 * t159 * t135 * t221 + t300) * t166;
t2 = (t41 * t147 - t296) * t148 - t150 * (t147 * t25 + t212);
t1 = (-t41 * t149 - t297) * t148 + (-0.2e1 * t147 * t122 * t223 + t301) * t150;
t186 = [t40 * t303 + t39 * t307 + t38 * t311 + t37 * t315, (-t28 * t245 - t27 * t246 - t26 * t247 - t25 * t248) * t184, t15 * t315 + t22 * t311 + t23 * t307 + t24 * t303 + (-t33 * t245 - t32 * t246 - t31 * t247 - t29 * t248) * t184, t14 * t315 + t19 * t311 + t20 * t307 + t21 * t303 + (-t36 * t245 - t35 * t246 - t34 * t247 - t30 * t248) * t184, (-t13 * t248 - t16 * t247 - t17 * t246 - t18 * t245) * t184, (-t10 * t247 - t11 * t246 - t12 * t245 - t9 * t248) * t318, (-t245 * t76 - t246 * t75 - t247 * t74 - t248 * t73) * t184, (-t245 * t72 - t246 * t71 - t247 * t70 - t248 * t69) * t184, 0, t1 * t315 + t3 * t311 + t4 * t307 + t5 * t303 + (-t111 * t317 - t112 * t313 - t113 * t309 - t114 * t305) * t184, t2 * t315 + t6 * t311 + t7 * t307 + t8 * t303 + (t111 * t224 + t112 * t222 + t113 * t220 + t114 * t218) * t184, 0, t92, -t91, -t115 * t77 + t116 * t78; t40 * t304 + t39 * t308 + t38 * t312 + t37 * t316, (t28 * t249 + t25 * t252 + t27 * t250 + t26 * t251) * t184, t15 * t316 + t22 * t312 + t23 * t308 + t24 * t304 + (t33 * t249 + t32 * t250 + t31 * t251 + t29 * t252) * t184, t14 * t316 + t19 * t312 + t20 * t308 + t21 * t304 + (t36 * t249 + t35 * t250 + t34 * t251 + t30 * t252) * t184, (t13 * t252 + t16 * t251 + t17 * t250 + t18 * t249) * t184, (t10 * t251 + t11 * t250 + t12 * t249 + t9 * t252) * t318, (t249 * t76 + t250 * t75 + t251 * t74 + t252 * t73) * t184, (t249 * t72 + t250 * t71 + t251 * t70 + t252 * t69) * t184, 0, t1 * t316 + t3 * t312 + t4 * t308 + t5 * t304 + (t107 * t317 + t108 * t313 + t109 * t309 + t110 * t305) * t184, t2 * t316 + t6 * t312 + t7 * t308 + t8 * t304 + (-t107 * t224 - t108 * t222 - t109 * t220 - t110 * t218) * t184, 0, t91, t92, t115 * t78 + t116 * t77; t40 * t234 + t39 * t237 + t38 * t240 + t37 * t244, (-t28 * t213 - t27 * t214 - t26 * t215 - t25 * t216) * t184, t15 * t244 + t22 * t240 + t23 * t237 + t24 * t234 + (-t33 * t213 - t32 * t214 - t31 * t215 - t29 * t216) * t184, t14 * t244 + t19 * t240 + t20 * t237 + t21 * t234 + (-t36 * t213 - t35 * t214 - t34 * t215 - t30 * t216) * t184, ((-t147 * t259 - t159 * t257 - t161 * t255 - t163 * t253) * t185 - t13 * t216 - t16 * t215 - t17 * t214 - t18 * t213) * t184, (-0.2e1 * t10 * t215 - 0.2e1 * t11 * t214 - 0.2e1 * t9 * t216 - 0.2e1 * t12 * t213 + (t141 * (-0.2e1 * t142 + t144) * t302 + t137 * (-0.2e1 * t138 + t140) * t306 + t133 * (-0.2e1 * t134 + t136) * t310 + t120 * (-0.2e1 * t121 + t123) * t314) * t185) * t184, ((t141 * t28 - t233 * t76) * t163 + (t137 * t27 - t236 * t75) * t161 + (t133 * t26 - t239 * t74) * t159 + (t120 * t25 - t243 * t73) * t147) * t184, (-t213 * t72 - t214 * t71 - t215 * t70 - t216 * t69 + t25 + t26 + t27 + t28) * t184, (t120 * t79 + t133 * t80 + t137 * t81 + t141 * t82) * t184, t1 * t244 + t3 * t240 + t4 * t237 + t5 * t234 + ((-g(3) * t169 + (-t279 * t33 + t36) * t163) * t141 + (-g(3) * t167 + (-t281 * t32 + t35) * t161) * t137 + (-g(3) * t165 + (-t283 * t31 + t34) * t159) * t133 + (-g(3) * t149 + (-t288 * t29 + t30) * t147) * t120) * t184, t2 * t244 + t6 * t240 + t7 * t237 + t8 * t234 + (t130 * t33 * t233 + t141 * (g(3) * t163 + t36 * t169) + t127 * t32 * t236 + t137 * (g(3) * t161 + t35 * t167) + t124 * t31 * t239 + t133 * (g(3) * t159 + t34 * t165) + t117 * t29 * t243 + t120 * (g(3) * t147 + t30 * t149)) * t184, 0, 0, 0, t156 - g(3); t46 * t37 + t51 * t38 + t52 * t39 + t53 * t40, (t25 * t260 + t28 * t254 + t27 * t256 + t26 * t258) * t184, t46 * t15 + t51 * t22 + t52 * t23 + t53 * t24 + (t33 * t254 + t32 * t256 + t31 * t258 + t29 * t260) * t184, t46 * t14 + t51 * t19 + t52 * t20 + t53 * t21 + (t36 * t254 + t35 * t256 + t34 * t258 + t30 * t260) * t184, (t13 * t260 + t16 * t258 + t17 * t256 + t18 * t254) * t184, (t10 * t258 + t11 * t256 + t12 * t254 + t9 * t260) * t318, (t254 * t76 + t256 * t75 + t258 * t74 + t260 * t73) * t184, (t254 * t72 + t256 * t71 + t258 * t70 + t260 * t69) * t184, 0, t46 * t1 + t51 * t3 + t52 * t4 + t53 * t5 + (t60 * t305 + t59 * t309 + t58 * t313 + t57 * t317) * t184, t46 * t2 + t51 * t6 + t52 * t7 + t53 * t8 + (-t218 * t60 - t220 * t59 - t222 * t58 - t224 * t57) * t184, t155, t77, -t78, 0;];
tauX_reg  = t186;
