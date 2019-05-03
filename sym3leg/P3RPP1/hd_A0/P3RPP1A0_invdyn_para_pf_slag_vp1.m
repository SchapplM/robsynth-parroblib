% Calculate vector of inverse dynamics forces for parallel robot
% P3RPP1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPP1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:51:56
% EndTime: 2019-05-03 14:52:00
% DurationCPUTime: 4.15s
% Computational Cost: add. (27509->391), mult. (37442->577), div. (2049->3), fcn. (14852->14), ass. (0->262)
t180 = pkin(1) + qJ(3,3);
t198 = xP(3);
t161 = sin(t198);
t162 = cos(t198);
t207 = koppelP(3,2);
t210 = koppelP(3,1);
t117 = t161 * t210 + t162 * t207;
t194 = xDP(3);
t196 = xDP(1);
t100 = -t117 * t194 + t196;
t201 = (qJ(3,3) ^ 2);
t218 = (pkin(1) ^ 2);
t271 = 1 + t218;
t273 = 2 * qJ(3,3);
t147 = pkin(1) * t273 + t201 + t271;
t202 = (qJ(2,3) ^ 2);
t286 = t147 + t202;
t144 = 1 / t286;
t186 = sin(qJ(1,3));
t189 = cos(qJ(1,3));
t248 = t180 * t189;
t241 = qJ(2,3) * t248;
t103 = t147 * t186 - t241;
t249 = t180 * t186;
t150 = qJ(2,3) * t249;
t106 = t147 * t189 + t150;
t177 = legFrame(3,3);
t153 = sin(t177);
t156 = cos(t177);
t70 = t103 * t156 + t106 * t153;
t73 = -t103 * t153 + t106 * t156;
t120 = -t161 * t207 + t162 * t210;
t195 = xDP(2);
t97 = t120 * t194 + t195;
t34 = (t100 * t73 + t70 * t97) * t144;
t295 = t180 * t34;
t181 = pkin(1) + qJ(3,2);
t208 = koppelP(2,2);
t211 = koppelP(2,1);
t118 = t161 * t211 + t162 * t208;
t101 = -t118 * t194 + t196;
t203 = qJ(3,2) ^ 2;
t275 = 2 * qJ(3,2);
t148 = pkin(1) * t275 + t203 + t271;
t204 = (qJ(2,2) ^ 2);
t285 = t148 + t204;
t145 = 1 / t285;
t187 = sin(qJ(1,2));
t190 = cos(qJ(1,2));
t246 = t181 * t190;
t242 = qJ(2,2) * t246;
t104 = t148 * t187 - t242;
t247 = t181 * t187;
t151 = qJ(2,2) * t247;
t107 = t148 * t190 + t151;
t178 = legFrame(2,3);
t154 = sin(t178);
t157 = cos(t178);
t71 = t104 * t157 + t107 * t154;
t74 = -t104 * t154 + t107 * t157;
t121 = -t161 * t208 + t162 * t211;
t98 = t121 * t194 + t195;
t35 = (t101 * t74 + t71 * t98) * t145;
t294 = t181 * t35;
t182 = pkin(1) + qJ(3,1);
t209 = koppelP(1,2);
t212 = koppelP(1,1);
t119 = t161 * t212 + t162 * t209;
t102 = -t119 * t194 + t196;
t205 = qJ(3,1) ^ 2;
t277 = 2 * qJ(3,1);
t149 = pkin(1) * t277 + t205 + t271;
t206 = (qJ(2,1) ^ 2);
t284 = t149 + t206;
t146 = 1 / t284;
t188 = sin(qJ(1,1));
t191 = cos(qJ(1,1));
t244 = t182 * t191;
t243 = qJ(2,1) * t244;
t105 = t149 * t188 - t243;
t245 = t182 * t188;
t152 = qJ(2,1) * t245;
t108 = t149 * t191 + t152;
t179 = legFrame(1,3);
t155 = sin(t179);
t158 = cos(t179);
t72 = t105 * t158 + t108 * t155;
t75 = -t105 * t155 + t108 * t158;
t122 = -t161 * t209 + t162 * t212;
t99 = t122 * t194 + t195;
t36 = (t102 * t75 + t72 * t99) * t146;
t293 = t182 * t36;
t166 = 1 + t202;
t167 = 1 + t204;
t168 = 1 + t206;
t137 = qJ(2,1) * t191 - t245;
t140 = qJ(2,1) * t188 + t244;
t93 = t137 * t158 - t140 * t155;
t96 = t137 * t155 + t140 * t158;
t289 = (t102 * t93 + t96 * t99) * t146;
t136 = qJ(2,2) * t190 - t247;
t139 = qJ(2,2) * t187 + t246;
t92 = t136 * t157 - t139 * t154;
t95 = t136 * t154 + t139 * t157;
t288 = (t101 * t92 + t95 * t98) * t145;
t135 = qJ(2,3) * t189 - t249;
t138 = qJ(2,3) * t186 + t248;
t91 = t135 * t156 - t138 * t153;
t94 = t135 * t153 + t138 * t156;
t287 = (t100 * t91 + t94 * t97) * t144;
t283 = -t201 - t202;
t282 = -t203 - t204;
t281 = -t205 - t206;
t280 = 2 * m(3);
t279 = 2 * pkin(1);
t278 = 0.2e1 * qJ(2,1);
t276 = 0.2e1 * qJ(2,2);
t274 = 0.2e1 * qJ(2,3);
t272 = -3 * t218;
t270 = (-qJ(2,3) - rSges(2,3)) * m(2);
t269 = (-qJ(2,2) - rSges(2,3)) * m(2);
t268 = (-qJ(2,1) - rSges(2,3)) * m(2);
t267 = m(3) * t144;
t266 = m(3) * t145;
t265 = m(3) * t146;
t172 = -qJ(2,3) - rSges(3,2);
t264 = m(3) * t172;
t174 = -qJ(2,2) - rSges(3,2);
t263 = m(3) * t174;
t176 = -qJ(2,1) - rSges(3,2);
t262 = m(3) * t176;
t109 = t166 * t186 + t241;
t112 = -t166 * t189 + t150;
t79 = t109 * t156 - t112 * t153;
t64 = t79 * t144 * t100;
t82 = t109 * t153 + t112 * t156;
t67 = t82 * t144 * t97;
t40 = t64 + t67;
t110 = t167 * t187 + t242;
t113 = -t167 * t190 + t151;
t80 = t110 * t157 - t113 * t154;
t65 = t80 * t145 * t101;
t83 = t110 * t154 + t113 * t157;
t68 = t83 * t145 * t98;
t41 = t65 + t68;
t111 = t168 * t188 + t243;
t114 = -t168 * t191 + t152;
t81 = t111 * t158 - t114 * t155;
t66 = t81 * t146 * t102;
t84 = t111 * t155 + t114 * t158;
t69 = t84 * t146 * t99;
t42 = t66 + t69;
t260 = qJ(2,1) * t42;
t258 = qJ(2,2) * t41;
t256 = qJ(2,3) * t40;
t255 = t144 * t287;
t254 = t145 * t288;
t253 = t146 * t289;
t252 = rSges(3,3) + qJ(3,1);
t251 = rSges(3,3) + qJ(3,2);
t250 = rSges(3,3) + qJ(3,3);
t240 = rSges(3,2) ^ 2 + (rSges(3,3) ^ 2) + t218;
t160 = (pkin(1) - rSges(2,2)) * m(2);
t163 = -pkin(1) - t250;
t141 = (m(3) * t163) - t160;
t164 = -pkin(1) - t251;
t142 = (m(3) * t164) - t160;
t165 = -pkin(1) - t252;
t143 = (m(3) * t165) - t160;
t239 = rSges(2,3) ^ 2 + t218 + ((-2 * pkin(1) + rSges(2,2)) * rSges(2,2));
t238 = -t264 - t270;
t237 = -t263 - t269;
t236 = -t262 - t268;
t229 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,1) + Icges(1,3);
t217 = pkin(1) * t218;
t200 = rSges(4,1);
t199 = rSges(4,2);
t197 = m(2) + m(3);
t193 = m(1) * rSges(1,1);
t192 = m(1) * rSges(1,2);
t185 = xDDP(1);
t184 = xDDP(2);
t183 = xDDP(3);
t170 = t194 ^ 2;
t134 = t192 - t236;
t133 = t192 - t237;
t132 = t192 - t238;
t131 = g(1) * t158 + g(2) * t155;
t130 = g(1) * t157 + g(2) * t154;
t129 = g(1) * t156 + g(2) * t153;
t128 = -g(1) * t155 + g(2) * t158;
t127 = -g(1) * t154 + g(2) * t157;
t126 = -g(1) * t153 + g(2) * t156;
t125 = -t143 + t193;
t124 = -t142 + t193;
t123 = -t141 + t193;
t116 = -t161 * t199 + t162 * t200;
t115 = t161 * t200 + t162 * t199;
t90 = -t119 * t183 - t122 * t170 + t185;
t89 = -t118 * t183 - t121 * t170 + t185;
t88 = -t117 * t183 - t120 * t170 + t185;
t87 = -t119 * t170 + t122 * t183 + t184;
t86 = -t118 * t170 + t121 * t183 + t184;
t85 = -t117 * t170 + t120 * t183 + t184;
t78 = (rSges(3,2) * t278 + (rSges(3,3) * t277) + (t252 * t279) + t240 - t281) * m(3) + (rSges(2,3) * t278 + t206 + t239) * m(2) + t229;
t77 = (rSges(3,2) * t276 + (rSges(3,3) * t275) + (t251 * t279) + t240 - t282) * m(3) + (rSges(2,3) * t276 + t204 + t239) * m(2) + t229;
t76 = (rSges(3,2) * t274 + (rSges(3,3) * t273) + (t250 * t279) + t240 - t283) * m(3) + (rSges(2,3) * t274 + t202 + t239) * m(2) + t229;
t63 = (t143 * t96 + t197 * t84) * t146;
t62 = (t142 * t95 + t197 * t83) * t145;
t61 = (t141 * t94 + t197 * t82) * t144;
t60 = (t143 * t93 + t197 * t81) * t146;
t59 = (t142 * t92 + t197 * t80) * t145;
t58 = (t141 * t91 + t197 * t79) * t144;
t57 = (-t176 * t96 + t72) * t265;
t56 = (-t174 * t95 + t71) * t266;
t55 = (-t172 * t94 + t70) * t267;
t54 = (-t176 * t93 + t75) * t265;
t53 = (-t174 * t92 + t74) * t266;
t52 = (-t172 * t91 + t73) * t267;
t51 = (-t119 * t93 + t122 * t96) * t146;
t50 = (-t118 * t92 + t121 * t95) * t145;
t49 = (-t117 * t91 + t120 * t94) * t144;
t45 = (-t119 * t81 + t122 * t84) * t146;
t44 = (-t118 * t80 + t121 * t83) * t145;
t43 = (-t117 * t79 + t120 * t82) * t144;
t39 = (-t119 * t75 + t122 * t72) * t146;
t38 = (-t118 * t74 + t121 * t71) * t145;
t37 = (-t117 * t73 + t120 * t70) * t144;
t33 = (t143 * t84 - t72 * t262 + t78 * t96) * t146;
t32 = (t142 * t83 - t71 * t263 + t77 * t95) * t145;
t31 = (t141 * t82 - t70 * t264 + t76 * t94) * t144;
t30 = (t143 * t81 - t75 * t262 + t78 * t93) * t146;
t29 = (t142 * t80 - t74 * t263 + t77 * t92) * t145;
t28 = (t141 * t79 - t73 * t264 + t76 * t91) * t144;
t27 = t143 * t51 + t197 * t45;
t26 = t142 * t50 + t197 * t44;
t25 = t141 * t49 + t197 * t43;
t24 = (-t176 * t51 + t39) * m(3);
t23 = (-t174 * t50 + t38) * m(3);
t22 = (-t172 * t49 + t37) * m(3);
t21 = t143 * t45 - t39 * t262 + t51 * t78;
t20 = t142 * t44 - t38 * t263 + t50 * t77;
t19 = t141 * t43 - t37 * t264 + t49 * t76;
t18 = 0.2e1 * (t260 + t293) * t253;
t17 = 0.2e1 * (t258 + t294) * t254;
t16 = 0.2e1 * (t256 + t295) * t255;
t15 = (0.2e1 * t182 * t260 - t36 + (-t206 - t168) * t36 - t284 * t289 * qJ(2,1)) * t253;
t14 = (0.2e1 * t181 * t258 - t35 + (-t204 - t167) * t35 - t285 * t288 * qJ(2,2)) * t254;
t13 = (0.2e1 * t180 * t256 - t34 + (-t202 - t166) * t34 - t286 * t287 * qJ(2,3)) * t255;
t12 = (0.2e1 * t149 * t42 - 0.2e1 * qJ(2,1) * t293 + (-t217 + (t272 + t281) * qJ(3,1) - qJ(3,1) + (-3 * t205 - t168) * pkin(1)) * t289) * t253;
t11 = (0.2e1 * t148 * t41 - 0.2e1 * qJ(2,2) * t294 + (-t217 + (t272 + t282) * qJ(3,2) - qJ(3,2) + (-3 * t203 - t167) * pkin(1)) * t288) * t254;
t10 = (0.2e1 * t147 * t40 - 0.2e1 * qJ(2,3) * t295 + (-t217 + (t272 + t283) * qJ(3,3) - qJ(3,3) + (-3 * t201 - t166) * pkin(1)) * t287) * t255;
t9 = -t143 * t18 - t289 * (t236 * t289 + t36 * t280) + (t128 * t191 - t131 * t188 - t15) * t197;
t8 = -t142 * t17 - t288 * (t237 * t288 + t35 * t280) + (t127 * t190 - t130 * t187 - t14) * t197;
t7 = -t141 * t16 - t287 * (t238 * t287 + t34 * t280) + (t126 * t189 - t129 * t186 - t13) * t197;
t6 = (t176 * t18 - t12 - t289 * (-t165 * t289 - 0.2e1 * t66 - 0.2e1 * t69) - t128 * t188 - t131 * t191) * m(3);
t5 = (t174 * t17 - t11 - t288 * (-t164 * t288 - 0.2e1 * t65 - 0.2e1 * t68) - t127 * t187 - t130 * t190) * m(3);
t4 = (t172 * t16 - t10 - t287 * (-t163 * t287 - 0.2e1 * t64 - 0.2e1 * t67) - t126 * t186 - t129 * t189) * m(3);
t3 = -t78 * t18 - t143 * t15 + t12 * t262 + 0.2e1 * t289 * ((-t165 * t36 - t176 * t42) * m(3) - t42 * t268) + (-t125 * t128 + t131 * t134) * t191 + (t125 * t131 + t128 * t134) * t188;
t2 = -t77 * t17 - t142 * t14 + t11 * t263 + 0.2e1 * t288 * ((-t164 * t35 - t174 * t41) * m(3) - t41 * t269) + (-t124 * t127 + t130 * t133) * t190 + (t124 * t130 + t127 * t133) * t187;
t1 = -t76 * t16 - t141 * t13 + t10 * t264 + 0.2e1 * t287 * ((-t163 * t34 - t172 * t40) * m(3) - t40 * t270) + (-t123 * t126 + t129 * t132) * t189 + (t123 * t129 + t126 * t132) * t186;
t46 = [(-t115 * t183 - t170 * t116 - g(1) + t185) * m(4) + ((t30 * t93 + t54 * t75 + t60 * t81) * t90 + (t30 * t96 + t54 * t72 + t60 * t84) * t87 + t93 * t3 + t81 * t9 + t75 * t6) * t146 + ((t29 * t92 + t53 * t74 + t59 * t80) * t89 + (t29 * t95 + t53 * t71 + t59 * t83) * t86 + t92 * t2 + t80 * t8 + t74 * t5) * t145 + ((t28 * t91 + t52 * t73 + t58 * t79) * t88 + (t28 * t94 + t52 * t70 + t58 * t82) * t85 + t91 * t1 + t79 * t7 + t73 * t4) * t144; (-t170 * t115 + t116 * t183 - g(2) + t184) * m(4) + ((t33 * t93 + t57 * t75 + t63 * t81) * t90 + (t33 * t96 + t57 * t72 + t63 * t84) * t87 + t96 * t3 + t84 * t9 + t72 * t6) * t146 + ((t32 * t92 + t56 * t74 + t62 * t80) * t89 + (t32 * t95 + t56 * t71 + t62 * t83) * t86 + t95 * t2 + t83 * t8 + t71 * t5) * t145 + ((t31 * t91 + t55 * t73 + t61 * t79) * t88 + (t31 * t94 + t55 * t70 + t61 * t82) * t85 + t94 * t1 + t82 * t7 + t70 * t4) * t144; Icges(4,3) * t183 + t49 * t1 + t50 * t2 + t51 * t3 + t37 * t4 + t38 * t5 + t39 * t6 + t43 * t7 + t44 * t8 + t45 * t9 + (-t115 * t185 + t116 * t184 + (t199 ^ 2 + t200 ^ 2) * t183 + (g(1) * t200 + g(2) * t199) * t161 + (g(1) * t199 - g(2) * t200) * t162) * m(4) + ((t21 * t93 + t24 * t75 + t27 * t81) * t90 + (t21 * t96 + t24 * t72 + t27 * t84) * t87) * t146 + ((t20 * t92 + t23 * t74 + t26 * t80) * t89 + (t20 * t95 + t23 * t71 + t26 * t83) * t86) * t145 + ((t19 * t91 + t22 * t73 + t25 * t79) * t88 + (t19 * t94 + t22 * t70 + t25 * t82) * t85) * t144;];
tauX  = t46;
