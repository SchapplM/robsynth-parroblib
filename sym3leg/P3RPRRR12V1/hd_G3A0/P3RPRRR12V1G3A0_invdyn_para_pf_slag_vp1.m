% Calculate vector of inverse dynamics forces for parallel robot
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:27:55
% EndTime: 2020-08-06 18:27:59
% DurationCPUTime: 3.86s
% Computational Cost: add. (9234->419), mult. (13509->756), div. (3243->7), fcn. (10260->18), ass. (0->256)
t148 = sin(qJ(3,3));
t112 = pkin(3) * t148 + qJ(2,3);
t109 = 0.1e1 / t112;
t163 = xDP(3);
t164 = xDP(2);
t165 = xDP(1);
t131 = 0.1e1 / t148;
t233 = t109 * t131;
t142 = legFrame(3,2);
t119 = sin(t142);
t122 = cos(t142);
t154 = cos(qJ(3,3));
t134 = t154 ^ 2;
t149 = sin(qJ(1,3));
t275 = (-pkin(5) - pkin(6));
t129 = pkin(1) - t275;
t155 = cos(qJ(1,3));
t181 = qJ(2,3) * t149 + t129 * t155;
t212 = t148 * t154;
t55 = t181 * t122 * t148 + t119 * t154 * qJ(2,3) + (t119 * t212 + (-t134 + 0.1e1) * t122 * t149) * pkin(3);
t224 = t122 * t154;
t58 = (pkin(3) * t224 - t181 * t119) * t148 + t149 * pkin(3) * (t154 - 0.1e1) * (t154 + 0.1e1) * t119 + qJ(2,3) * t224;
t82 = t112 * t155 - t129 * t149;
t22 = t82 * t109 * t163 + (t164 * t58 + t165 * t55) * t233;
t286 = 0.2e1 * t22;
t150 = sin(qJ(3,2));
t113 = pkin(3) * t150 + qJ(2,2);
t110 = 0.1e1 / t113;
t132 = 0.1e1 / t150;
t232 = t110 * t132;
t143 = legFrame(2,2);
t120 = sin(t143);
t123 = cos(t143);
t156 = cos(qJ(3,2));
t135 = t156 ^ 2;
t151 = sin(qJ(1,2));
t157 = cos(qJ(1,2));
t182 = qJ(2,2) * t151 + t129 * t157;
t211 = t150 * t156;
t56 = t182 * t123 * t150 + t120 * t156 * qJ(2,2) + (t120 * t211 + (-t135 + 0.1e1) * t123 * t151) * pkin(3);
t222 = t123 * t156;
t59 = (pkin(3) * t222 - t182 * t120) * t150 + t151 * pkin(3) * (t156 - 0.1e1) * (t156 + 0.1e1) * t120 + qJ(2,2) * t222;
t83 = t113 * t157 - t129 * t151;
t23 = t83 * t110 * t163 + (t164 * t59 + t165 * t56) * t232;
t285 = 0.2e1 * t23;
t152 = sin(qJ(3,1));
t114 = pkin(3) * t152 + qJ(2,1);
t111 = 0.1e1 / t114;
t133 = 0.1e1 / t152;
t231 = t111 * t133;
t144 = legFrame(1,2);
t121 = sin(t144);
t124 = cos(t144);
t158 = cos(qJ(3,1));
t136 = t158 ^ 2;
t153 = sin(qJ(1,1));
t159 = cos(qJ(1,1));
t183 = qJ(2,1) * t153 + t129 * t159;
t210 = t152 * t158;
t57 = t183 * t124 * t152 + t121 * t158 * qJ(2,1) + (t121 * t210 + (-t136 + 0.1e1) * t124 * t153) * pkin(3);
t220 = t124 * t158;
t60 = (pkin(3) * t220 - t183 * t121) * t152 + t153 * pkin(3) * (t158 - 0.1e1) * (t158 + 0.1e1) * t121 + qJ(2,1) * t220;
t84 = t114 * t159 - t129 * t153;
t24 = t84 * t111 * t163 + (t164 * t60 + t165 * t57) * t231;
t284 = 0.2e1 * t24;
t67 = (-t149 * t163 + (-t119 * t164 + t122 * t165) * t155) * t109;
t283 = 0.2e1 * t67;
t68 = (-t151 * t163 + (-t120 * t164 + t123 * t165) * t157) * t110;
t282 = 0.2e1 * t68;
t69 = (-t153 * t163 + (-t121 * t164 + t124 * t165) * t159) * t111;
t281 = 0.2e1 * t69;
t280 = 2 * pkin(1);
t279 = -0.2e1 * pkin(3);
t278 = 2 * rSges(2,3);
t171 = (rSges(3,2) ^ 2);
t173 = (rSges(3,1) ^ 2);
t101 = (m(3) * (-t171 + t173) - Icges(3,1) + Icges(3,2));
t277 = 2 * t101;
t274 = m(3) * rSges(3,2);
t117 = (rSges(3,1) * t274 - Icges(3,4));
t276 = 2 * t117;
t98 = rSges(3,1) * t154 - rSges(3,2) * t148;
t273 = m(3) * t98;
t99 = rSges(3,1) * t156 - rSges(3,2) * t150;
t272 = m(3) * t99;
t177 = 0.1e1 / pkin(3);
t217 = t131 * t177;
t79 = (-t119 * t165 - t122 * t164) * t217;
t271 = pkin(3) * t79;
t215 = t132 * t177;
t80 = (-t120 * t165 - t123 * t164) * t215;
t270 = pkin(3) * t80;
t213 = t133 * t177;
t81 = (-t121 * t165 - t124 * t164) * t213;
t269 = pkin(3) * t81;
t160 = (rSges(3,3) + pkin(5));
t100 = rSges(3,1) * t158 - rSges(3,2) * t152;
t268 = m(3) * t100;
t267 = m(3) * (-pkin(1) - t160);
t266 = t67 * t79;
t265 = t68 * t80;
t264 = t69 * t81;
t263 = rSges(3,1) * t148;
t262 = rSges(3,1) * t150;
t261 = rSges(3,1) * t152;
t260 = qJ(2,1) * t81;
t259 = qJ(2,2) * t80;
t258 = qJ(2,3) * t79;
t257 = t129 * t67;
t256 = t129 * t68;
t255 = t129 * t69;
t103 = -rSges(3,2) * t267 - Icges(3,6);
t104 = -rSges(3,1) * t267 - Icges(3,5);
t73 = t103 * t148 - t104 * t154;
t254 = (-t149 * t73 + t82 * t273) * t233;
t253 = t131 * t55;
t252 = t131 * t58;
t166 = m(2) + m(3);
t97 = -(pkin(1) - rSges(2,2)) * m(2) + t267;
t61 = (-t149 * t97 + t166 * t82) * t109;
t251 = t131 * t61;
t250 = t131 * t97;
t249 = t131 * t98;
t74 = t103 * t150 - t104 * t156;
t248 = (-t151 * t74 + t83 * t272) * t232;
t247 = t132 * t56;
t246 = t132 * t59;
t62 = (-t151 * t97 + t166 * t83) * t110;
t245 = t132 * t62;
t244 = t132 * t97;
t243 = t132 * t99;
t75 = t103 * t152 - t104 * t158;
t242 = (-t153 * t75 + t84 * t268) * t231;
t241 = t133 * t57;
t240 = t133 * t60;
t63 = (-t153 * t97 + t166 * t84) * t111;
t239 = t133 * t63;
t238 = t133 * t97;
t237 = qJ(2,1) * t152;
t236 = qJ(2,2) * t150;
t235 = qJ(2,3) * t148;
t234 = t100 * t133;
t230 = t117 * t134;
t229 = t117 * t135;
t228 = t117 * t136;
t227 = t119 * t155;
t226 = t120 * t157;
t225 = t121 * t159;
t223 = t122 * t155;
t221 = t123 * t157;
t219 = t124 * t159;
t218 = t131 * t166;
t216 = t132 * t166;
t214 = t133 * t166;
t201 = -t274 / 0.2e1;
t209 = rSges(3,1) * t201 + Icges(3,4) / 0.2e1;
t208 = m(3) * t249;
t207 = m(3) * t243;
t206 = t154 * t271;
t205 = t156 * t270;
t204 = t158 * t269;
t203 = m(3) * t234;
t202 = m(3) * rSges(3,1) / 0.2e1;
t200 = t119 * t217;
t199 = t120 * t215;
t198 = t121 * t213;
t197 = t122 * t217;
t196 = t123 * t215;
t195 = t124 * t213;
t194 = t177 * t208;
t193 = t177 * t207;
t192 = t177 * t203;
t178 = pkin(1) ^ 2;
t191 = t173 + t178 + (t280 + t160) * t160;
t190 = rSges(2,3) ^ 2 + t178 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t107 = qJ(2,1) * m(3) + m(2) * (rSges(2,3) + qJ(2,1));
t106 = qJ(2,2) * m(3) + m(2) * (rSges(2,3) + qJ(2,2));
t105 = qJ(2,3) * m(3) + m(2) * (rSges(2,3) + qJ(2,3));
t90 = g(1) * t122 - g(2) * t119;
t189 = g(3) * t155 + t149 * t90;
t91 = g(1) * t123 - g(2) * t120;
t188 = g(3) * t157 + t151 * t91;
t92 = g(1) * t124 - g(2) * t121;
t187 = g(3) * t159 + t153 * t92;
t186 = rSges(3,2) * t154 + t263;
t185 = rSges(3,2) * t156 + t262;
t184 = rSges(3,2) * t158 + t261;
t180 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,2) + Icges(1,3);
t176 = pkin(3) ^ 2;
t179 = -(pkin(6) ^ 2) + (t275 * t280) - t176 - t178 + ((-2 * pkin(6) - pkin(5)) * pkin(5));
t169 = qJ(2,1) ^ 2;
t168 = qJ(2,2) ^ 2;
t167 = qJ(2,3) ^ 2;
t162 = m(1) * rSges(1,2);
t147 = xDDP(1);
t146 = xDDP(2);
t145 = xDDP(3);
t108 = (t171 + t173) * m(3) + Icges(3,3);
t96 = -t107 + t162;
t95 = -t106 + t162;
t94 = -t105 + t162;
t93 = (t173 / 0.2e1 - t171 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t89 = g(1) * t121 + g(2) * t124;
t88 = g(1) * t120 + g(2) * t123;
t87 = g(1) * t119 + g(2) * t122;
t86 = m(1) * rSges(1,1) - t97;
t85 = g(3) * t86;
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t66 = t69 ^ 2;
t65 = t68 ^ 2;
t64 = t67 ^ 2;
t45 = t66 + t78;
t44 = t65 + t77;
t43 = t64 + t76;
t42 = -t101 * t136 + t210 * t276 + (qJ(2,1) * t278 + t169 + t190) * m(2) + (0.2e1 * qJ(2,1) * t184 + t169 + t191) * m(3) + t180;
t41 = -t101 * t135 + t211 * t276 + (qJ(2,2) * t278 + t168 + t190) * m(2) + (0.2e1 * qJ(2,2) * t185 + t168 + t191) * m(3) + t180;
t40 = -t101 * t134 + t212 * t276 + (qJ(2,3) * t278 + t167 + t190) * m(2) + (0.2e1 * qJ(2,3) * t186 + t167 + t191) * m(3) + t180;
t39 = -t124 * t192 + (t60 * t214 - t97 * t225) * t111;
t38 = -t123 * t193 + (t59 * t216 - t97 * t226) * t110;
t37 = -t122 * t194 + (t58 * t218 - t97 * t227) * t109;
t36 = -t121 * t192 + (t57 * t214 + t97 * t219) * t111;
t35 = -t120 * t193 + (t56 * t216 + t97 * t221) * t110;
t34 = -t119 * t194 + (t55 * t218 + t97 * t223) * t109;
t33 = (-t153 * t42 + t84 * t97) * t111;
t32 = (-t151 * t41 + t83 * t97) * t110;
t31 = (-t149 * t40 + t82 * t97) * t109;
t30 = -t108 * t195 + (t60 * t203 - t75 * t225) * t111;
t29 = -t108 * t196 + (t59 * t207 - t74 * t226) * t110;
t28 = -t108 * t197 + (t58 * t208 - t73 * t227) * t109;
t27 = -t108 * t198 + (t57 * t203 + t75 * t219) * t111;
t26 = -t108 * t199 + (t56 * t207 + t74 * t221) * t110;
t25 = -t108 * t200 + (t55 * t208 + t73 * t223) * t109;
t21 = -t75 * t195 + (-t42 * t225 + t60 * t238) * t111;
t20 = -t74 * t196 + (-t41 * t226 + t59 * t244) * t110;
t19 = -t73 * t197 + (-t40 * t227 + t58 * t250) * t109;
t18 = -t75 * t198 + (t42 * t219 + t57 * t238) * t111;
t17 = -t74 * t199 + (t41 * t221 + t56 * t244) * t110;
t16 = -t73 * t200 + (t40 * t223 + t55 * t250) * t109;
t15 = (0.2e1 * t204 + t284 - t255) * t69 * t111;
t14 = (0.2e1 * t205 + t285 - t256) * t68 * t110;
t13 = (0.2e1 * t206 + t286 - t257) * t67 * t109;
t12 = (((t158 * t255 - t269) * t152 - t260) * t133 * t269 + t255 * t284 + (t129 * t204 + (t136 * t176 + t237 * t279 - t169 + t179) * t69) * t69) * t111;
t11 = (((t156 * t256 - t270) * t150 - t259) * t132 * t270 + t256 * t285 + (t129 * t205 + (t135 * t176 + t236 * t279 - t168 + t179) * t68) * t68) * t110;
t10 = (((t154 * t257 - t271) * t148 - t258) * t131 * t271 + t257 * t286 + (t129 * t206 + (t134 * t176 + t235 * t279 - t167 + t179) * t67) * t67) * t109;
t9 = -t75 * t15 - t12 * t268 - t108 * t158 * t78 * t133 - 0.2e1 * t66 * (t228 + (qJ(2,1) * t202 + t93 * t152) * t158 + t201 * t237 + t209) - ((t187 * rSges(3,1) - rSges(3,2) * t89) * t158 - t152 * (rSges(3,1) * t89 + t187 * rSges(3,2))) * m(3);
t8 = -t74 * t14 - t11 * t272 - t108 * t156 * t77 * t132 - 0.2e1 * t65 * (t229 + (qJ(2,2) * t202 + t93 * t150) * t156 + t201 * t236 + t209) - ((t188 * rSges(3,1) - rSges(3,2) * t88) * t156 - t150 * (rSges(3,1) * t88 + t188 * rSges(3,2))) * m(3);
t7 = -t73 * t13 - t10 * t273 - t108 * t154 * t76 * t131 - 0.2e1 * t64 * (t230 + (qJ(2,3) * t202 + t93 * t148) * t154 + t201 * t235 + t209) - ((t189 * rSges(3,1) - rSges(3,2) * t87) * t154 - t148 * (rSges(3,1) * t87 + t189 * rSges(3,2))) * m(3);
t6 = -t107 * t66 - t97 * t15 + (-t12 - t187) * t166 + (-t45 * t261 + (-rSges(3,2) * t45 - t78 * t234) * t158) * m(3);
t5 = -t106 * t65 - t97 * t14 + (-t11 - t188) * t166 + (-t44 * t262 + (-rSges(3,2) * t44 - t77 * t243) * t156) * m(3);
t4 = -t105 * t64 - t97 * t13 + (-t10 - t189) * t166 + (-t43 * t263 + (-rSges(3,2) * t43 - t76 * t249) * t154) * m(3);
t3 = -t42 * t15 - t97 * t12 + 0.4e1 * t228 * t264 + t104 * t78 * t152 + (t107 * t24 - t117 * t81) * t281 + (t92 * t96 + t85) * t159 + (-g(3) * t96 + t86 * t92) * t153 + (t152 * t264 * t277 + (-t133 * t75 + t103) * t78) * t158 + (((rSges(3,1) * t260 + rSges(3,2) * t24) * t158 + (rSges(3,1) * t24 - rSges(3,2) * t260) * t152) * t281 + (t153 * g(3) - t159 * t92) * t184) * m(3);
t2 = -t41 * t14 - t97 * t11 + 0.4e1 * t229 * t265 + t104 * t77 * t150 + (t106 * t23 - t117 * t80) * t282 + (t91 * t95 + t85) * t157 + (-g(3) * t95 + t86 * t91) * t151 + (t150 * t265 * t277 + (-t132 * t74 + t103) * t77) * t156 + (((rSges(3,1) * t259 + rSges(3,2) * t23) * t156 + (rSges(3,1) * t23 - rSges(3,2) * t259) * t150) * t282 + (t151 * g(3) - t157 * t91) * t185) * m(3);
t1 = -t40 * t13 - t97 * t10 + 0.4e1 * t230 * t266 + t104 * t76 * t148 + (t105 * t22 - t117 * t79) * t283 + (t90 * t94 + t85) * t155 + (-g(3) * t94 + t86 * t90) * t149 + (t148 * t266 * t277 + (-t131 * t73 + t103) * t76) * t154 + (((rSges(3,1) * t258 + rSges(3,2) * t22) * t154 + (rSges(3,1) * t22 - rSges(3,2) * t258) * t148) * t283 + (t149 * g(3) - t155 * t90) * t186) * m(3);
t46 = [(-g(1) + t147) * m(4) + ((t18 * t219 + t36 * t241) * t147 + (-t18 * t225 + t36 * t240) * t146 + (-t153 * t18 + t36 * t84) * t145 + t3 * t219 + t6 * t241) * t111 + ((-t124 * t146 * t27 + (-t147 * t27 - t9) * t121) * t133 + (-t123 * t146 * t26 + (-t147 * t26 - t8) * t120) * t132 + (-t122 * t146 * t25 + (-t147 * t25 - t7) * t119) * t131) * t177 + ((t17 * t221 + t35 * t247) * t147 + (-t17 * t226 + t35 * t246) * t146 + (-t151 * t17 + t35 * t83) * t145 + t2 * t221 + t5 * t247) * t110 + ((t16 * t223 + t34 * t253) * t147 + (-t16 * t227 + t34 * t252) * t146 + (-t149 * t16 + t34 * t82) * t145 + t1 * t223 + t4 * t253) * t109; (-g(2) + t146) * m(4) + ((t21 * t219 + t39 * t241) * t147 + (-t21 * t225 + t39 * t240) * t146 + (-t153 * t21 + t39 * t84) * t145 - t3 * t225 + t6 * t240) * t111 + ((-t121 * t147 * t30 + (-t146 * t30 - t9) * t124) * t133 + (-t120 * t147 * t29 + (-t146 * t29 - t8) * t123) * t132 + (-t119 * t147 * t28 + (-t146 * t28 - t7) * t122) * t131) * t177 + ((t20 * t221 + t38 * t247) * t147 + (-t20 * t226 + t38 * t246) * t146 + (-t151 * t20 + t38 * t83) * t145 - t2 * t226 + t5 * t246) * t110 + ((t19 * t223 + t37 * t253) * t147 + (-t19 * t227 + t37 * t252) * t146 + (-t149 * t19 + t37 * t82) * t145 - t1 * t227 + t4 * t252) * t109; (-g(3) + t145) * m(4) + ((-t119 * t254 - t120 * t248 - t121 * t242) * t147 + (-t122 * t254 - t123 * t248 - t124 * t242) * t146) * t177 + ((t33 * t219 + t57 * t239) * t147 + (-t33 * t225 + t60 * t239) * t146 + (-t153 * t33 + t63 * t84) * t145 - t153 * t3 + t84 * t6) * t111 + ((t32 * t221 + t56 * t245) * t147 + (-t32 * t226 + t59 * t245) * t146 + (-t151 * t32 + t62 * t83) * t145 - t151 * t2 + t83 * t5) * t110 + ((t31 * t223 + t55 * t251) * t147 + (-t31 * t227 + t58 * t251) * t146 + (-t149 * t31 + t61 * t82) * t145 - t149 * t1 + t82 * t4) * t109;];
tauX  = t46;
