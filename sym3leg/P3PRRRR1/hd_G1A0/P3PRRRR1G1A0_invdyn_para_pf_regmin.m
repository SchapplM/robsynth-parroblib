% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR1G1P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR1G1P3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1P3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:09
% EndTime: 2020-03-09 20:34:10
% DurationCPUTime: 1.83s
% Computational Cost: add. (2612->248), mult. (6836->511), div. (2670->17), fcn. (7164->18), ass. (0->226)
t146 = sin(qJ(2,3));
t119 = 0.1e1 / t146 ^ 2;
t139 = legFrame(3,3);
t112 = sin(t139);
t115 = cos(t139);
t157 = xDP(2);
t158 = xDP(1);
t100 = t112 * t158 - t115 * t157;
t151 = cos(qJ(3,3));
t145 = sin(qJ(3,3));
t152 = cos(qJ(2,3));
t229 = t145 * t152;
t82 = (-t112 * t157 - t115 * t158) * t151 - t100 * t229;
t79 = t82 ^ 2;
t275 = t119 * t79;
t148 = sin(qJ(2,2));
t122 = 0.1e1 / t148 ^ 2;
t140 = legFrame(2,3);
t113 = sin(t140);
t116 = cos(t140);
t101 = t113 * t158 - t116 * t157;
t153 = cos(qJ(3,2));
t147 = sin(qJ(3,2));
t154 = cos(qJ(2,2));
t228 = t147 * t154;
t83 = (-t113 * t157 - t116 * t158) * t153 - t101 * t228;
t80 = t83 ^ 2;
t274 = t122 * t80;
t150 = sin(qJ(2,1));
t125 = 0.1e1 / t150 ^ 2;
t141 = legFrame(1,3);
t114 = sin(t141);
t117 = cos(t141);
t102 = t114 * t158 - t117 * t157;
t155 = cos(qJ(3,1));
t149 = sin(qJ(3,1));
t156 = cos(qJ(2,1));
t227 = t149 * t156;
t84 = (-t114 * t157 - t117 * t158) * t155 - t102 * t227;
t81 = t84 ^ 2;
t273 = t125 * t81;
t159 = 0.1e1 / pkin(2);
t118 = 0.1e1 / t146;
t121 = 0.1e1 / t148;
t124 = 0.1e1 / t150;
t127 = 0.1e1 / t151;
t131 = 0.1e1 / t153;
t135 = 0.1e1 / t155;
t160 = 0.1e1 / pkin(2) ^ 2;
t170 = t151 ^ 2;
t128 = 0.1e1 / t170;
t173 = t153 ^ 2;
t132 = 0.1e1 / t173;
t176 = t155 ^ 2;
t136 = 0.1e1 / t176;
t142 = xDDP(3);
t272 = t142 - g(3);
t109 = t115 * g(1) + t112 * g(2);
t143 = xDDP(2);
t144 = xDDP(1);
t129 = t127 * t128;
t232 = t129 * t159;
t270 = t118 * t275;
t226 = t151 * t152;
t88 = t112 * t145 + t115 * t226;
t94 = t112 * t226 - t115 * t145;
t97 = t100 ^ 2;
t52 = t232 * t270 + (t97 * t232 + (t143 * t94 + t144 * t88) * t127) * t118 + t272;
t196 = t109 * t146 + t52 * t152;
t271 = t118 * t196;
t110 = t116 * g(1) + t113 * g(2);
t133 = t131 * t132;
t231 = t133 * t159;
t268 = t121 * t274;
t225 = t153 * t154;
t90 = t113 * t147 + t116 * t225;
t95 = t113 * t225 - t116 * t147;
t98 = t101 ^ 2;
t53 = t231 * t268 + (t98 * t231 + (t143 * t95 + t144 * t90) * t131) * t121 + t272;
t195 = t110 * t148 + t53 * t154;
t269 = t121 * t195;
t111 = t117 * g(1) + t114 * g(2);
t137 = t135 * t136;
t230 = t137 * t159;
t266 = t124 * t273;
t224 = t155 * t156;
t92 = t114 * t149 + t117 * t224;
t96 = t114 * t224 - t117 * t149;
t99 = t102 ^ 2;
t54 = t230 * t266 + (t99 * t230 + (t143 * t96 + t144 * t92) * t135) * t124 + t272;
t194 = t111 * t150 + t54 * t156;
t267 = t124 * t194;
t130 = 0.1e1 / t170 ^ 2;
t265 = t130 * t79;
t134 = 0.1e1 / t173 ^ 2;
t264 = t134 * t80;
t138 = 0.1e1 / t176 ^ 2;
t263 = t138 * t81;
t237 = t118 * t128;
t247 = t100 * t118;
t85 = -t112 * t229 - t115 * t151;
t89 = -t151 * t112 + t115 * t229;
t28 = ((t143 * t89 + t144 * t85) * t237 + (-(-t100 * t145 * t146 + t118 * t152 * t82) * t119 * t82 - (t100 * t152 - t145 * t82) * t247) * t127 * t232) * t159;
t262 = t145 * t28;
t235 = t121 * t132;
t246 = t101 * t121;
t86 = -t113 * t228 - t116 * t153;
t91 = -t153 * t113 + t116 * t228;
t29 = ((t143 * t91 + t144 * t86) * t235 + (-(-t101 * t147 * t148 + t121 * t154 * t83) * t122 * t83 - (t101 * t154 - t147 * t83) * t246) * t131 * t231) * t159;
t261 = t147 * t29;
t233 = t124 * t136;
t245 = t102 * t124;
t87 = -t114 * t227 - t117 * t155;
t93 = -t155 * t114 + t117 * t227;
t30 = ((t143 * t93 + t144 * t87) * t233 + (-(-t102 * t149 * t150 + t124 * t156 * t84) * t125 * t84 - (t102 * t156 - t149 * t84) * t245) * t135 * t230) * t159;
t260 = t149 * t30;
t259 = t160 * t97;
t258 = t160 * t98;
t257 = t160 * t99;
t256 = t28 * t146;
t255 = t29 * t148;
t254 = t30 * t150;
t202 = t145 * t259;
t76 = t129 * t202 + (t112 * t144 - t115 * t143) * t159 * t127;
t253 = t76 * t145;
t252 = t76 * t151;
t201 = t147 * t258;
t77 = t133 * t201 + (t113 * t144 - t116 * t143) * t159 * t131;
t251 = t77 * t147;
t250 = t77 * t153;
t200 = t149 * t257;
t78 = t137 * t200 + (t114 * t144 - t117 * t143) * t159 * t135;
t249 = t78 * t149;
t248 = t78 * t155;
t244 = t112 * t127;
t243 = t113 * t131;
t242 = t114 * t135;
t241 = t115 * t127;
t240 = t116 * t131;
t239 = t117 * t135;
t238 = t118 * t127;
t236 = t121 * t131;
t234 = t124 * t135;
t190 = t160 * t82 * t247;
t218 = t119 * t265;
t73 = (t128 * t97 + t218) * t160;
t1 = t28 * t226 + (-t73 * t151 - t253) * t146 - 0.2e1 * t129 * t190 * t229;
t189 = t160 * t83 * t246;
t212 = t122 * t264;
t74 = (t132 * t98 + t212) * t160;
t2 = t29 * t225 + (-t74 * t153 - t251) * t148 - 0.2e1 * t133 * t189 * t228;
t188 = t160 * t84 * t245;
t206 = t125 * t263;
t75 = (t136 * t99 + t206) * t160;
t3 = t30 * t224 + (-t75 * t155 - t249) * t150 - 0.2e1 * t137 * t188 * t227;
t223 = t88 * t238;
t222 = t94 * t238;
t221 = t85 * t237;
t220 = t89 * t237;
t219 = t118 * t265;
t217 = t90 * t236;
t216 = t95 * t236;
t215 = t86 * t235;
t214 = t91 * t235;
t213 = t121 * t264;
t211 = t92 * t234;
t210 = t96 * t234;
t209 = t87 * t233;
t208 = t93 * t233;
t207 = t124 * t263;
t205 = t127 * t262;
t204 = t131 * t261;
t203 = t135 * t260;
t199 = 0.2e1 * (t151 * t262 + (0.2e1 * t127 - t129) * t190) * t237;
t198 = 0.2e1 * (t153 * t261 + (0.2e1 * t131 - t133) * t189) * t235;
t197 = 0.2e1 * (t155 * t260 + (0.2e1 * t135 - t137) * t188) * t233;
t13 = 0.2e1 * t136 * t188 + t260;
t6 = -t156 * t13 + (t75 * t149 - t248) * t150;
t14 = 0.2e1 * t128 * t190 + t262;
t4 = -t152 * t14 + (t73 * t145 - t252) * t146;
t15 = 0.2e1 * t132 * t189 + t261;
t5 = -t154 * t15 + (t74 * t147 - t250) * t148;
t193 = t109 * t152 - t146 * t52;
t192 = t110 * t154 - t148 * t53;
t191 = t111 * t156 - t150 * t54;
t187 = t145 * t196 * t237;
t186 = t145 * t218;
t185 = t152 * t218;
t184 = t147 * t195 * t235;
t183 = t147 * t212;
t182 = t154 * t212;
t181 = t149 * t194 * t233;
t180 = t149 * t206;
t179 = t156 * t206;
t161 = t159 * t160;
t108 = t114 * g(1) - t117 * g(2);
t107 = t113 * g(1) - t116 * g(2);
t106 = t112 * g(1) - t115 * g(2);
t72 = t135 * t257 + t249;
t71 = t131 * t258 + t251;
t70 = t127 * t259 + t253;
t66 = -t136 * t200 + t248;
t65 = -t132 * t201 + t250;
t64 = -t128 * t202 + t252;
t63 = (-0.2e1 * t136 + t138) * t160 * t273;
t62 = (-0.2e1 * t132 + t134) * t160 * t274;
t61 = (-0.2e1 * t128 + t130) * t160 * t275;
t36 = t149 * t108 + t155 * t191;
t35 = -t155 * t108 + t149 * t191;
t34 = t147 * t107 + t153 * t192;
t33 = -t153 * t107 + t147 * t192;
t32 = t145 * t106 + t151 * t193;
t31 = -t151 * t106 + t145 * t193;
t27 = t156 * t30;
t26 = t154 * t29;
t25 = t152 * t28;
t21 = -t160 * t207 + t27;
t20 = -t160 * t213 + t26;
t19 = -t160 * t219 + t25;
t18 = -t160 * t179 - t254;
t17 = -t160 * t182 - t255;
t16 = -t160 * t185 - t256;
t12 = t13 * t149;
t11 = t15 * t147;
t10 = t14 * t145;
t7 = [t54 * t211 + t53 * t217 + t52 * t223, (t30 * t209 + t29 * t215 + t28 * t221) * t159, t19 * t223 + t20 * t217 + t21 * t211 + (t194 * t209 + t195 * t215 + t196 * t221) * t159, t16 * t223 + t17 * t217 + t18 * t211 + (t191 * t209 + t192 * t215 + t193 * t221) * t159, (-t112 * t186 - t113 * t183 - t114 * t180) * t161 + (t10 * t221 + t11 * t215 + t12 * t209) * t159, (t197 * t87 + t198 * t86 + t199 * t85 + t63 * t242 + t62 * t243 + t61 * t244) * t159, (t112 * t205 + t113 * t204 + t114 * t203 + t209 * t72 + t215 * t71 + t221 * t70) * t159, (t112 * t28 + t113 * t29 + t114 * t30 + t209 * t66 + t215 * t65 + t221 * t64) * t159, (t78 * t242 + t77 * t243 + t76 * t244) * t159, t1 * t223 + t2 * t217 + t3 * t211 + ((t114 * t35 + t87 * t267) * t135 + (t113 * t33 + t86 * t269) * t131 + (t112 * t31 + t85 * t271) * t127) * t159, t4 * t223 + t5 * t217 + t6 * t211 + (-t181 * t87 - t184 * t86 - t187 * t85 + t36 * t242 + t34 * t243 + t32 * t244) * t159, t144 - g(1); t54 * t210 + t53 * t216 + t52 * t222, (t30 * t208 + t29 * t214 + t28 * t220) * t159, t19 * t222 + t20 * t216 + t21 * t210 + (t194 * t208 + t195 * t214 + t196 * t220) * t159, t16 * t222 + t17 * t216 + t18 * t210 + (t191 * t208 + t192 * t214 + t193 * t220) * t159, (t115 * t186 + t116 * t183 + t117 * t180) * t161 + (t10 * t220 + t11 * t214 + t12 * t208) * t159, (t197 * t93 + t198 * t91 + t199 * t89 - t63 * t239 - t62 * t240 - t61 * t241) * t159, (-t115 * t205 - t116 * t204 - t117 * t203 + t208 * t72 + t214 * t71 + t220 * t70) * t159, (-t115 * t28 - t116 * t29 - t117 * t30 + t208 * t66 + t214 * t65 + t220 * t64) * t159, (-t78 * t239 - t77 * t240 - t76 * t241) * t159, t1 * t222 + t2 * t216 + t3 * t210 + ((-t117 * t35 + t93 * t267) * t135 + (-t116 * t33 + t91 * t269) * t131 + (-t115 * t31 + t89 * t271) * t127) * t159, t4 * t222 + t5 * t216 + t6 * t210 + (-t181 * t93 - t184 * t91 - t187 * t89 - t36 * t239 - t34 * t240 - t32 * t241) * t159, t143 - g(2); -(3 * g(3)) + (3 * t142) + (t211 + t217 + t223) * t144 + (t210 + t216 + t222) * t143 + ((t124 * t99 + t266) * t137 + (t121 * t98 + t268) * t133 + (t118 * t97 + t270) * t129) * t159, 0, t25 + t26 + t27 + (-t207 - t213 - t219) * t160, -t256 - t255 - t254 + (-t179 - t182 - t185) * t160, 0, 0, 0, 0, 0, t3 + t2 + t1, t4 + t5 + t6, t272;];
tauX_reg  = t7;
