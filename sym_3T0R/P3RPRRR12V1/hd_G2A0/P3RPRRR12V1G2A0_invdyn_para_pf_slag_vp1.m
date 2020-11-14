% Calculate vector of inverse dynamics forces for parallel robot
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:23
% EndTime: 2020-08-06 18:24:27
% DurationCPUTime: 3.50s
% Computational Cost: add. (9234->419), mult. (13350->754), div. (3243->7), fcn. (10101->18), ass. (0->260)
t151 = sin(qJ(3,3));
t115 = pkin(3) * t151 + qJ(2,3);
t112 = 0.1e1 / t115;
t166 = xDP(3);
t167 = xDP(2);
t168 = xDP(1);
t134 = 0.1e1 / t151;
t233 = t112 * t134;
t145 = legFrame(3,2);
t125 = cos(t145);
t157 = cos(qJ(3,3));
t122 = sin(t145);
t226 = t122 * t157;
t158 = cos(qJ(1,3));
t269 = pkin(3) * t158;
t278 = (-pkin(5) - pkin(6));
t132 = pkin(1) - t278;
t152 = sin(qJ(1,3));
t93 = qJ(2,3) * t158 - t132 * t152;
t55 = (pkin(3) * t226 - t125 * t93) * t151 + (t157 - 0.1e1) * (t157 + 0.1e1) * t125 * t269 + qJ(2,3) * t226;
t137 = t157 ^ 2;
t220 = t125 * t157;
t58 = (pkin(3) * t220 + t122 * t93) * t151 + (-t137 + 0.1e1) * t122 * t269 + qJ(2,3) * t220;
t82 = t115 * t152 + t132 * t158;
t22 = t82 * t112 * t166 + (t167 * t58 + t168 * t55) * t233;
t290 = 0.2e1 * t22;
t153 = sin(qJ(3,2));
t116 = pkin(3) * t153 + qJ(2,2);
t113 = 0.1e1 / t116;
t135 = 0.1e1 / t153;
t232 = t113 * t135;
t146 = legFrame(2,2);
t126 = cos(t146);
t159 = cos(qJ(3,2));
t123 = sin(t146);
t224 = t123 * t159;
t160 = cos(qJ(1,2));
t268 = pkin(3) * t160;
t154 = sin(qJ(1,2));
t94 = qJ(2,2) * t160 - t132 * t154;
t56 = (pkin(3) * t224 - t126 * t94) * t153 + (t159 - 0.1e1) * (t159 + 0.1e1) * t126 * t268 + qJ(2,2) * t224;
t138 = t159 ^ 2;
t218 = t126 * t159;
t59 = (pkin(3) * t218 + t123 * t94) * t153 + (-t138 + 0.1e1) * t123 * t268 + qJ(2,2) * t218;
t83 = t116 * t154 + t132 * t160;
t23 = t83 * t113 * t166 + (t167 * t59 + t168 * t56) * t232;
t289 = 0.2e1 * t23;
t155 = sin(qJ(3,1));
t117 = pkin(3) * t155 + qJ(2,1);
t114 = 0.1e1 / t117;
t136 = 0.1e1 / t155;
t231 = t114 * t136;
t147 = legFrame(1,2);
t127 = cos(t147);
t161 = cos(qJ(3,1));
t124 = sin(t147);
t222 = t124 * t161;
t162 = cos(qJ(1,1));
t267 = pkin(3) * t162;
t156 = sin(qJ(1,1));
t95 = qJ(2,1) * t162 - t132 * t156;
t57 = (pkin(3) * t222 - t127 * t95) * t155 + (t161 - 0.1e1) * (t161 + 0.1e1) * t127 * t267 + qJ(2,1) * t222;
t139 = t161 ^ 2;
t216 = t127 * t161;
t60 = (pkin(3) * t216 + t124 * t95) * t155 + (-t139 + 0.1e1) * t124 * t267 + qJ(2,1) * t216;
t84 = t117 * t156 + t132 * t162;
t24 = t84 * t114 * t166 + (t167 * t60 + t168 * t57) * t231;
t288 = 0.2e1 * t24;
t67 = (t158 * t166 + (-t122 * t167 + t125 * t168) * t152) * t112;
t287 = 0.2e1 * t67;
t68 = (t160 * t166 + (-t123 * t167 + t126 * t168) * t154) * t113;
t286 = 0.2e1 * t68;
t69 = (t162 * t166 + (-t124 * t167 + t127 * t168) * t156) * t114;
t285 = 0.2e1 * t69;
t284 = 2 * pkin(1);
t283 = -0.2e1 * pkin(3);
t282 = 2 * rSges(2,3);
t281 = 0.2e1 * t151;
t280 = 0.2e1 * t153;
t279 = 0.2e1 * t155;
t277 = m(3) * rSges(3,2);
t180 = 0.1e1 / pkin(3);
t214 = t134 * t180;
t79 = (-t122 * t168 - t125 * t167) * t214;
t276 = pkin(3) * t79;
t212 = t135 * t180;
t80 = (-t123 * t168 - t126 * t167) * t212;
t275 = pkin(3) * t80;
t210 = t136 * t180;
t81 = (-t124 * t168 - t127 * t167) * t210;
t274 = pkin(3) * t81;
t163 = (rSges(3,3) + pkin(5));
t101 = rSges(3,1) * t157 - rSges(3,2) * t151;
t273 = m(3) * t101;
t102 = rSges(3,1) * t159 - rSges(3,2) * t153;
t272 = m(3) * t102;
t103 = rSges(3,1) * t161 - rSges(3,2) * t155;
t271 = m(3) * t103;
t270 = m(3) * (-pkin(1) - t163);
t266 = t67 * t79;
t265 = t68 * t80;
t264 = t69 * t81;
t263 = rSges(3,1) * t151;
t262 = rSges(3,1) * t153;
t261 = rSges(3,1) * t155;
t260 = qJ(2,1) * t81;
t259 = qJ(2,2) * t80;
t258 = qJ(2,3) * t79;
t257 = t132 * t67;
t256 = t132 * t68;
t255 = t132 * t69;
t106 = -rSges(3,2) * t270 - Icges(3,6);
t107 = -rSges(3,1) * t270 - Icges(3,5);
t73 = t106 * t151 - t107 * t157;
t254 = (t158 * t73 + t82 * t273) * t233;
t253 = t134 * t55;
t252 = t134 * t58;
t100 = -((pkin(1) - rSges(2,2)) * m(2)) + t270;
t169 = m(2) + m(3);
t61 = (t100 * t158 + t169 * t82) * t112;
t251 = t134 * t61;
t74 = t106 * t153 - t107 * t159;
t250 = (t160 * t74 + t83 * t272) * t232;
t249 = t135 * t56;
t248 = t135 * t59;
t62 = (t100 * t160 + t169 * t83) * t113;
t247 = t135 * t62;
t75 = t106 * t155 - t107 * t161;
t246 = (t162 * t75 + t84 * t271) * t231;
t245 = t136 * t57;
t244 = t136 * t60;
t63 = (t100 * t162 + t169 * t84) * t114;
t243 = t136 * t63;
t242 = qJ(2,1) * t155;
t241 = qJ(2,2) * t153;
t240 = qJ(2,3) * t151;
t239 = t100 * t134;
t238 = t100 * t135;
t237 = t100 * t136;
t236 = t101 * t134;
t235 = t102 * t135;
t234 = t103 * t136;
t120 = rSges(3,1) * t277 - Icges(3,4);
t230 = t120 * t137;
t229 = t120 * t138;
t228 = t120 * t139;
t227 = t122 * t152;
t225 = t123 * t154;
t223 = t124 * t156;
t221 = t125 * t152;
t219 = t126 * t154;
t217 = t127 * t156;
t215 = t134 * t169;
t213 = t135 * t169;
t211 = t136 * t169;
t201 = -t277 / 0.2e1;
t209 = rSges(3,1) * t201 + Icges(3,4) / 0.2e1;
t208 = t157 * t276;
t207 = t159 * t275;
t206 = t161 * t274;
t205 = m(3) * t236;
t204 = m(3) * t235;
t203 = m(3) * t234;
t202 = m(3) * rSges(3,1) / 0.2e1;
t200 = t122 * t214;
t199 = t123 * t212;
t198 = t124 * t210;
t197 = t125 * t214;
t196 = t126 * t212;
t195 = t127 * t210;
t194 = t180 * t205;
t193 = t180 * t204;
t192 = t180 * t203;
t176 = rSges(3,1) ^ 2;
t181 = pkin(1) ^ 2;
t191 = t176 + t181 + ((t284 + t163) * t163);
t190 = rSges(2,3) ^ 2 + t181 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t174 = rSges(3,2) ^ 2;
t104 = m(3) * (-t174 + t176) - Icges(3,1) + Icges(3,2);
t110 = qJ(2,1) * m(3) + m(2) * (rSges(2,3) + qJ(2,1));
t109 = qJ(2,2) * m(3) + m(2) * (rSges(2,3) + qJ(2,2));
t108 = qJ(2,3) * m(3) + m(2) * (rSges(2,3) + qJ(2,3));
t90 = g(1) * t125 - g(2) * t122;
t189 = g(3) * t152 - t158 * t90;
t91 = g(1) * t126 - g(2) * t123;
t188 = g(3) * t154 - t160 * t91;
t92 = g(1) * t127 - g(2) * t124;
t187 = g(3) * t156 - t162 * t92;
t186 = rSges(3,2) * t157 + t263;
t185 = rSges(3,2) * t159 + t262;
t184 = rSges(3,2) * t161 + t261;
t183 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,2) + Icges(1,3);
t179 = pkin(3) ^ 2;
t182 = -(pkin(6) ^ 2) + (t278 * t284) - t179 - t181 + ((-2 * pkin(6) - pkin(5)) * pkin(5));
t172 = qJ(2,1) ^ 2;
t171 = qJ(2,2) ^ 2;
t170 = qJ(2,3) ^ 2;
t165 = m(1) * rSges(1,2);
t150 = xDDP(1);
t149 = xDDP(2);
t148 = xDDP(3);
t111 = (t174 + t176) * m(3) + Icges(3,3);
t99 = -t110 + t165;
t98 = -t109 + t165;
t97 = -t108 + t165;
t96 = (t176 / 0.2e1 - t174 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t89 = g(1) * t124 + g(2) * t127;
t88 = g(1) * t123 + g(2) * t126;
t87 = g(1) * t122 + g(2) * t125;
t86 = m(1) * rSges(1,1) - t100;
t85 = t86 * g(3);
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t66 = t69 ^ 2;
t65 = t68 ^ 2;
t64 = t67 ^ 2;
t45 = t66 + t78;
t44 = t65 + t77;
t43 = t64 + t76;
t42 = -t104 * t139 + t120 * t161 * t279 + (qJ(2,1) * t282 + t172 + t190) * m(2) + (0.2e1 * qJ(2,1) * t184 + t172 + t191) * m(3) + t183;
t41 = -t104 * t138 + t120 * t159 * t280 + (qJ(2,2) * t282 + t171 + t190) * m(2) + (0.2e1 * qJ(2,2) * t185 + t171 + t191) * m(3) + t183;
t40 = -t104 * t137 + t120 * t157 * t281 + (qJ(2,3) * t282 + t170 + t190) * m(2) + (0.2e1 * qJ(2,3) * t186 + t170 + t191) * m(3) + t183;
t39 = -t127 * t192 + (-t100 * t223 + t60 * t211) * t114;
t38 = -t126 * t193 + (-t100 * t225 + t59 * t213) * t113;
t37 = -t125 * t194 + (-t100 * t227 + t58 * t215) * t112;
t36 = -t124 * t192 + (t100 * t217 + t57 * t211) * t114;
t35 = -t123 * t193 + (t100 * t219 + t56 * t213) * t113;
t34 = -t122 * t194 + (t100 * t221 + t55 * t215) * t112;
t33 = (t100 * t84 + t162 * t42) * t114;
t32 = (t100 * t83 + t160 * t41) * t113;
t31 = (t100 * t82 + t158 * t40) * t112;
t30 = -t111 * t195 + (t60 * t203 - t75 * t223) * t114;
t29 = -t111 * t196 + (t59 * t204 - t74 * t225) * t113;
t28 = -t111 * t197 + (t58 * t205 - t73 * t227) * t112;
t27 = -t111 * t198 + (t57 * t203 + t75 * t217) * t114;
t26 = -t111 * t199 + (t56 * t204 + t74 * t219) * t113;
t25 = -t111 * t200 + (t55 * t205 + t73 * t221) * t112;
t21 = -t75 * t195 + (-t42 * t223 + t60 * t237) * t114;
t20 = -t74 * t196 + (-t41 * t225 + t59 * t238) * t113;
t19 = -t73 * t197 + (-t40 * t227 + t58 * t239) * t112;
t18 = -t75 * t198 + (t42 * t217 + t57 * t237) * t114;
t17 = -t74 * t199 + (t41 * t219 + t56 * t238) * t113;
t16 = -t73 * t200 + (t40 * t221 + t55 * t239) * t112;
t15 = (0.2e1 * t206 + t288 - t255) * t69 * t114;
t14 = (0.2e1 * t207 + t289 - t256) * t68 * t113;
t13 = (0.2e1 * t208 + t290 - t257) * t67 * t112;
t12 = (((t161 * t255 - t274) * t155 - t260) * t136 * t274 + t255 * t288 + (t132 * t206 + (t139 * t179 + t242 * t283 - t172 + t182) * t69) * t69) * t114;
t11 = (((t159 * t256 - t275) * t153 - t259) * t135 * t275 + t256 * t289 + (t132 * t207 + (t138 * t179 + t241 * t283 - t171 + t182) * t68) * t68) * t113;
t10 = (((t157 * t257 - t276) * t151 - t258) * t134 * t276 + t257 * t290 + (t132 * t208 + (t137 * t179 + t240 * t283 - t170 + t182) * t67) * t67) * t112;
t9 = -t75 * t15 - t12 * t271 - t111 * t161 * t78 * t136 - 0.2e1 * (t228 + (qJ(2,1) * t202 + t96 * t155) * t161 + t201 * t242 + t209) * t66 - ((t187 * rSges(3,1) - rSges(3,2) * t89) * t161 - t155 * (rSges(3,1) * t89 + t187 * rSges(3,2))) * m(3);
t8 = -t74 * t14 - t11 * t272 - t111 * t159 * t77 * t135 - 0.2e1 * (t229 + (qJ(2,2) * t202 + t96 * t153) * t159 + t201 * t241 + t209) * t65 - ((t188 * rSges(3,1) - rSges(3,2) * t88) * t159 - t153 * (rSges(3,1) * t88 + t188 * rSges(3,2))) * m(3);
t7 = -t73 * t13 - t10 * t273 - t111 * t157 * t76 * t134 - 0.2e1 * (t230 + (qJ(2,3) * t202 + t96 * t151) * t157 + t201 * t240 + t209) * t64 - ((t189 * rSges(3,1) - rSges(3,2) * t87) * t157 - t151 * (rSges(3,1) * t87 + t189 * rSges(3,2))) * m(3);
t6 = -t100 * t15 - t66 * t110 + (-t12 - t187) * t169 + (-t45 * t261 + (-rSges(3,2) * t45 - t78 * t234) * t161) * m(3);
t5 = -t100 * t14 - t65 * t109 + (-t11 - t188) * t169 + (-t44 * t262 + (-rSges(3,2) * t44 - t77 * t235) * t159) * m(3);
t4 = -t100 * t13 - t64 * t108 + (-t10 - t189) * t169 + (-t43 * t263 + (-rSges(3,2) * t43 - t76 * t236) * t157) * m(3);
t3 = -t42 * t15 - t100 * t12 + 0.4e1 * t228 * t264 + t107 * t78 * t155 + (t110 * t24 - t120 * t81) * t285 + (g(3) * t99 - t86 * t92) * t162 + (t92 * t99 + t85) * t156 + (t104 * t264 * t279 + (-t136 * t75 + t106) * t78) * t161 + (((rSges(3,1) * t260 + rSges(3,2) * t24) * t161 + (rSges(3,1) * t24 - rSges(3,2) * t260) * t155) * t285 + (-t162 * g(3) - t156 * t92) * t184) * m(3);
t2 = -t41 * t14 - t100 * t11 + 0.4e1 * t229 * t265 + t107 * t77 * t153 + (t109 * t23 - t120 * t80) * t286 + (g(3) * t98 - t86 * t91) * t160 + (t91 * t98 + t85) * t154 + (t104 * t265 * t280 + (-t135 * t74 + t106) * t77) * t159 + (((rSges(3,1) * t259 + rSges(3,2) * t23) * t159 + (rSges(3,1) * t23 - rSges(3,2) * t259) * t153) * t286 + (-t160 * g(3) - t154 * t91) * t185) * m(3);
t1 = -t40 * t13 - t100 * t10 + 0.4e1 * t230 * t266 + t107 * t76 * t151 + (t108 * t22 - t120 * t79) * t287 + (g(3) * t97 - t86 * t90) * t158 + (t90 * t97 + t85) * t152 + (t104 * t266 * t281 + (-t134 * t73 + t106) * t76) * t157 + (((rSges(3,1) * t258 + rSges(3,2) * t22) * t157 + (rSges(3,1) * t22 - rSges(3,2) * t258) * t151) * t287 + (-t158 * g(3) - t152 * t90) * t186) * m(3);
t46 = [(-g(1) + t150) * m(4) + ((t18 * t217 + t36 * t245) * t150 + (-t18 * t223 + t36 * t244) * t149 + (t162 * t18 + t36 * t84) * t148 + t3 * t217 + t6 * t245) * t114 + ((-t127 * t149 * t27 + (-t150 * t27 - t9) * t124) * t136 + (-t126 * t149 * t26 + (-t150 * t26 - t8) * t123) * t135 + (-t125 * t149 * t25 + (-t150 * t25 - t7) * t122) * t134) * t180 + ((t17 * t219 + t35 * t249) * t150 + (-t17 * t225 + t35 * t248) * t149 + (t160 * t17 + t35 * t83) * t148 + t2 * t219 + t5 * t249) * t113 + ((t16 * t221 + t34 * t253) * t150 + (-t16 * t227 + t34 * t252) * t149 + (t158 * t16 + t34 * t82) * t148 + t1 * t221 + t4 * t253) * t112; (-g(2) + t149) * m(4) + ((t21 * t217 + t39 * t245) * t150 + (-t21 * t223 + t39 * t244) * t149 + (t162 * t21 + t39 * t84) * t148 - t3 * t223 + t6 * t244) * t114 + ((-t124 * t150 * t30 + (-t149 * t30 - t9) * t127) * t136 + (-t123 * t150 * t29 + (-t149 * t29 - t8) * t126) * t135 + (-t122 * t150 * t28 + (-t149 * t28 - t7) * t125) * t134) * t180 + ((t20 * t219 + t38 * t249) * t150 + (-t20 * t225 + t38 * t248) * t149 + (t160 * t20 + t38 * t83) * t148 - t2 * t225 + t5 * t248) * t113 + ((t19 * t221 + t37 * t253) * t150 + (-t19 * t227 + t37 * t252) * t149 + (t158 * t19 + t37 * t82) * t148 - t1 * t227 + t4 * t252) * t112; (-g(3) + t148) * m(4) + ((-t122 * t254 - t123 * t250 - t124 * t246) * t150 + (-t125 * t254 - t126 * t250 - t127 * t246) * t149) * t180 + ((t33 * t217 + t57 * t243) * t150 + (-t33 * t223 + t60 * t243) * t149 + (t162 * t33 + t63 * t84) * t148 + t162 * t3 + t84 * t6) * t114 + ((t32 * t219 + t56 * t247) * t150 + (-t32 * t225 + t59 * t247) * t149 + (t160 * t32 + t62 * t83) * t148 + t160 * t2 + t83 * t5) * t113 + ((t31 * t221 + t55 * t251) * t150 + (-t31 * t227 + t58 * t251) * t149 + (t158 * t31 + t61 * t82) * t148 + t158 * t1 + t82 * t4) * t112;];
tauX  = t46;
