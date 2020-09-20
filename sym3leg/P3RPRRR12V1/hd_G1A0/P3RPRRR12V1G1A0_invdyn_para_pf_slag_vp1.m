% Calculate vector of inverse dynamics forces for parallel robot
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
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:16
% EndTime: 2020-08-06 18:21:19
% DurationCPUTime: 2.76s
% Computational Cost: add. (6150->360), mult. (7818->624), div. (1476->11), fcn. (6312->18), ass. (0->225)
t256 = 2 * qJ(2,1);
t255 = 2 * qJ(2,2);
t254 = 2 * qJ(2,3);
t165 = (rSges(3,2) ^ 2);
t167 = (rSges(3,1) ^ 2);
t88 = (m(3) * (-t165 + t167) - Icges(3,1) + Icges(3,2));
t257 = -2 * t88;
t141 = sin(qJ(3,3));
t118 = 0.1e1 / t141;
t143 = sin(qJ(3,2));
t121 = 0.1e1 / t143;
t145 = sin(qJ(3,1));
t124 = 0.1e1 / t145;
t119 = 0.1e1 / t141 ^ 2;
t122 = 0.1e1 / t143 ^ 2;
t125 = 0.1e1 / t145 ^ 2;
t252 = m(3) * rSges(3,2);
t103 = (rSges(3,1) * t252 - Icges(3,4));
t253 = 2 * t103;
t147 = cos(qJ(3,3));
t85 = rSges(3,1) * t147 - rSges(3,2) * t141;
t251 = m(3) * t85;
t149 = cos(qJ(3,2));
t86 = rSges(3,1) * t149 - rSges(3,2) * t143;
t250 = m(3) * t86;
t151 = cos(qJ(3,1));
t87 = rSges(3,1) * t151 - rSges(3,2) * t145;
t249 = m(3) * t87;
t248 = rSges(3,2) * g(3);
t157 = xDP(3);
t204 = t147 * t118;
t196 = t157 * t204;
t158 = xDP(2);
t159 = xDP(1);
t135 = legFrame(3,3);
t105 = sin(t135);
t108 = cos(t135);
t116 = pkin(1) + pkin(5) + pkin(6);
t142 = sin(qJ(1,3));
t148 = cos(qJ(1,3));
t99 = pkin(3) * t141 + qJ(2,3);
t61 = t116 * t148 + t142 * t99;
t64 = t116 * t142 - t148 * t99;
t49 = t105 * t61 + t108 * t64;
t96 = 0.1e1 / t99;
t236 = t49 * t96;
t46 = -t105 * t64 + t108 * t61;
t239 = t46 * t96;
t233 = t158 * t236 + t159 * t239;
t19 = t196 + t233;
t67 = -t105 * t142 + t108 * t148;
t70 = t105 * t148 + t108 * t142;
t43 = (t158 * t70 + t159 * t67) * t96;
t227 = t116 * t43;
t242 = t43 * t96;
t10 = (t19 - t196 - t227 + t233) * t242;
t120 = t118 * t119;
t127 = t147 ^ 2;
t230 = rSges(3,1) * t141;
t185 = rSges(3,2) * t147 + t230;
t169 = 0.1e1 / pkin(3);
t201 = t157 * t169;
t188 = qJ(2,3) * t118 * t201;
t205 = t157 ^ 2 / pkin(3) ^ 2;
t161 = qJ(2,3) ^ 2;
t179 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,2) + Icges(1,3);
t171 = (pkin(1) ^ 2);
t186 = rSges(2,3) ^ 2 + t171 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t153 = (rSges(3,3) + pkin(5));
t187 = t167 + t171 + (2 * pkin(1) + t153) * t153;
t28 = -t88 * t127 + t141 * t147 * t253 + ((rSges(2,3) * t254 + t161 + t186) * m(2)) + (t185 * t254 + t161 + t187) * m(3) + t179;
t244 = m(3) * (-pkin(1) - t153);
t90 = -rSges(3,2) * t244 - Icges(3,6);
t91 = -rSges(3,1) * t244 - Icges(3,5);
t58 = t141 * t90 - t147 * t91;
t168 = pkin(3) ^ 2;
t178 = pkin(3) * t116 * t201;
t199 = -t116 ^ 2 - t168;
t214 = t19 * t116;
t7 = -((t118 * t157 + t147 * t227) * t141 + t188) * t96 * t119 * t157 + (((-t178 * t204 + t214) * t141 + ((t127 * t168 - t161 + t199) * t141 + (t127 - 0.1e1) * t254 * pkin(3)) * t43) * t118 + t214) * t242;
t84 = -((pkin(1) - rSges(2,2)) * m(2)) + t244;
t73 = m(1) * rSges(1,1) - t84;
t74 = -g(1) * t105 + g(2) * t108;
t77 = g(1) * t108 + g(2) * t105;
t156 = m(1) * rSges(1,2);
t92 = qJ(2,3) * m(3) + (m(2) * (rSges(2,3) + qJ(2,3)));
t81 = t156 - t92;
t247 = t96 * (-t28 * t10 - t84 * t7 + (-t73 * t74 + t77 * t81) * t148 + (t73 * t77 + t74 * t81) * t142 + (0.2e1 * t19 * t92 + (t147 * t257 + (-0.4e1 * t127 + 0.2e1) * t118 * t103) * t201) * t43 + (-t120 * t147 * t58 + (t141 * t91 + t147 * t90) * t119) * t205 + (0.2e1 * ((-rSges(3,1) * t188 + rSges(3,2) * t19) * t147 + (rSges(3,1) * t19 + rSges(3,2) * t188) * t141) * t43 + (-t142 * t74 - t148 * t77) * t185) * m(3));
t203 = t149 * t121;
t194 = t157 * t203;
t136 = legFrame(2,3);
t106 = sin(t136);
t109 = cos(t136);
t100 = pkin(3) * t143 + qJ(2,2);
t144 = sin(qJ(1,2));
t150 = cos(qJ(1,2));
t62 = t100 * t144 + t116 * t150;
t65 = -t100 * t150 + t116 * t144;
t50 = t106 * t62 + t109 * t65;
t97 = 0.1e1 / t100;
t235 = t50 * t97;
t47 = -t106 * t65 + t109 * t62;
t238 = t47 * t97;
t232 = t158 * t235 + t159 * t238;
t20 = t194 + t232;
t68 = -t106 * t144 + t109 * t150;
t71 = t106 * t150 + t109 * t144;
t44 = (t158 * t71 + t159 * t68) * t97;
t226 = t116 * t44;
t241 = t44 * t97;
t11 = (-t194 + t20 - t226 + t232) * t241;
t123 = t121 * t122;
t128 = t149 ^ 2;
t229 = rSges(3,1) * t143;
t184 = rSges(3,2) * t149 + t229;
t189 = qJ(2,2) * t121 * t201;
t162 = qJ(2,2) ^ 2;
t29 = -t88 * t128 + t143 * t149 * t253 + ((rSges(2,3) * t255 + t162 + t186) * m(2)) + (t184 * t255 + t162 + t187) * m(3) + t179;
t59 = t143 * t90 - t149 * t91;
t75 = -g(1) * t106 + g(2) * t109;
t78 = g(1) * t109 + g(2) * t106;
t213 = t20 * t116;
t8 = -((t121 * t157 + t149 * t226) * t143 + t189) * t97 * t122 * t157 + (((-t178 * t203 + t213) * t143 + ((t128 * t168 - t162 + t199) * t143 + (t128 - 0.1e1) * t255 * pkin(3)) * t44) * t121 + t213) * t241;
t93 = qJ(2,2) * m(3) + (m(2) * (rSges(2,3) + qJ(2,2)));
t82 = t156 - t93;
t246 = t97 * (-t29 * t11 - t84 * t8 + (-t73 * t75 + t78 * t82) * t150 + (t73 * t78 + t75 * t82) * t144 + (0.2e1 * t20 * t93 + (t149 * t257 + (-0.4e1 * t128 + 0.2e1) * t121 * t103) * t201) * t44 + (-t123 * t149 * t59 + (t143 * t91 + t149 * t90) * t122) * t205 + (0.2e1 * ((-rSges(3,1) * t189 + rSges(3,2) * t20) * t149 + (rSges(3,1) * t20 + rSges(3,2) * t189) * t143) * t44 + (-t144 * t75 - t150 * t78) * t184) * m(3));
t202 = t151 * t124;
t192 = t157 * t202;
t137 = legFrame(1,3);
t107 = sin(t137);
t110 = cos(t137);
t101 = pkin(3) * t145 + qJ(2,1);
t146 = sin(qJ(1,1));
t152 = cos(qJ(1,1));
t63 = t101 * t146 + t116 * t152;
t66 = -t101 * t152 + t116 * t146;
t51 = t107 * t63 + t110 * t66;
t98 = 0.1e1 / t101;
t234 = t51 * t98;
t48 = -t107 * t66 + t110 * t63;
t237 = t48 * t98;
t231 = t158 * t234 + t159 * t237;
t21 = t192 + t231;
t69 = -t107 * t146 + t110 * t152;
t72 = t107 * t152 + t110 * t146;
t45 = (t158 * t72 + t159 * t69) * t98;
t225 = t116 * t45;
t240 = t45 * t98;
t12 = (-t192 + t21 - t225 + t231) * t240;
t126 = t124 * t125;
t129 = t151 ^ 2;
t228 = rSges(3,1) * t145;
t183 = rSges(3,2) * t151 + t228;
t190 = qJ(2,1) * t124 * t201;
t163 = qJ(2,1) ^ 2;
t30 = -t88 * t129 + t145 * t151 * t253 + ((rSges(2,3) * t256 + t163 + t186) * m(2)) + (t183 * t256 + t163 + t187) * m(3) + t179;
t60 = t145 * t90 - t151 * t91;
t76 = -g(1) * t107 + g(2) * t110;
t79 = g(1) * t110 + g(2) * t107;
t94 = qJ(2,1) * m(3) + (m(2) * (rSges(2,3) + qJ(2,1)));
t83 = t156 - t94;
t212 = t21 * t116;
t9 = -((t124 * t157 + t151 * t225) * t145 + t190) * t98 * t125 * t157 + (((-t178 * t202 + t212) * t145 + ((t129 * t168 - t163 + t199) * t145 + (t129 - 0.1e1) * t256 * pkin(3)) * t45) * t124 + t212) * t240;
t245 = t98 * (-t30 * t12 - t84 * t9 + (-t73 * t76 + t79 * t83) * t152 + (t73 * t79 + t76 * t83) * t146 + (0.2e1 * t21 * t94 + (t151 * t257 + (-0.4e1 * t129 + 0.2e1) * t124 * t103) * t201) * t45 + (-t126 * t151 * t60 + (t145 * t91 + t151 * t90) * t125) * t205 + (0.2e1 * ((-rSges(3,1) * t190 + rSges(3,2) * t21) * t151 + (rSges(3,1) * t21 + rSges(3,2) * t190) * t145) * t45 + (-t146 * t76 - t152 * t79) * t183) * m(3));
t243 = m(3) * t169;
t139 = xDDP(2);
t224 = t139 * t96;
t223 = t139 * t97;
t222 = t139 * t98;
t140 = xDDP(1);
t221 = t140 * t96;
t220 = t140 * t97;
t219 = t140 * t98;
t95 = (t165 + t167) * m(3) + Icges(3,3);
t218 = t169 * t95;
t217 = t169 * t96;
t216 = t169 * t97;
t215 = t169 * t98;
t138 = xDDP(3);
t211 = t118 * t138;
t210 = t118 * t169;
t209 = t121 * t138;
t208 = t121 * t169;
t207 = t124 * t138;
t206 = t124 * t169;
t197 = -t252 / 0.2e1;
t200 = rSges(3,1) * t197 + Icges(3,4) / 0.2e1;
t198 = m(3) * rSges(3,1) / 0.2e1;
t195 = t120 * t205;
t193 = t123 * t205;
t191 = t126 * t205;
t182 = t142 * t77 - t148 * t74;
t181 = t144 * t78 - t150 * t75;
t180 = t146 * t79 - t152 * t76;
t160 = m(2) + m(3);
t155 = rSges(3,1) * g(3);
t80 = (t167 / 0.2e1 - t165 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t57 = (t151 * t160 - t87 * t243) * t124;
t56 = (t149 * t160 - t86 * t243) * t121;
t55 = (t147 * t160 - t85 * t243) * t118;
t54 = (t151 * t84 - t169 * t60) * t124;
t53 = (t149 * t84 - t169 * t59) * t121;
t52 = (t147 * t84 - t169 * t58) * t118;
t42 = t45 ^ 2;
t41 = t44 ^ 2;
t40 = t43 ^ 2;
t33 = t125 * t205 + t42;
t32 = t122 * t205 + t41;
t31 = t119 * t205 + t40;
t27 = (t160 * t51 + t72 * t84) * t98;
t26 = (t160 * t50 + t71 * t84) * t97;
t25 = (t160 * t49 + t70 * t84) * t96;
t24 = (t160 * t48 + t69 * t84) * t98;
t23 = (t160 * t47 + t68 * t84) * t97;
t22 = (t160 * t46 + t67 * t84) * t96;
t18 = (t30 * t72 + t51 * t84) * t98;
t17 = (t29 * t71 + t50 * t84) * t97;
t16 = (t28 * t70 + t49 * t84) * t96;
t15 = (t30 * t69 + t48 * t84) * t98;
t14 = (t29 * t68 + t47 * t84) * t97;
t13 = (t28 * t67 + t46 * t84) * t96;
t6 = -t84 * t12 - t94 * t42 + (-t180 - t9) * t160 + (-t33 * t228 + (-rSges(3,2) * t33 - t87 * t191) * t151) * m(3);
t5 = -t84 * t11 - t93 * t41 + (-t181 - t8) * t160 + (-t32 * t229 + (-rSges(3,2) * t32 - t86 * t193) * t149) * m(3);
t4 = -t84 * t10 - t92 * t40 + (-t182 - t7) * t160 + (-t31 * t230 + (-rSges(3,2) * t31 - t85 * t195) * t147) * m(3);
t1 = [(t15 * t69 + t24 * t48) * t219 + (t15 * t72 + t24 * t51) * t222 + (t24 * t151 - (t48 * t249 + t60 * t69) * t215) * t207 + t69 * t245 + t6 * t237 + (t14 * t68 + t23 * t47) * t220 + (t14 * t71 + t23 * t50) * t223 + (t23 * t149 - (t47 * t250 + t59 * t68) * t216) * t209 + t68 * t246 + t5 * t238 + (t13 * t67 + t22 * t46) * t221 + (t13 * t70 + t22 * t49) * t224 + (t22 * t147 - (t46 * t251 + t58 * t67) * t217) * t211 + t67 * t247 + t4 * t239 + (-g(1) + t140) * m(4); (t18 * t69 + t27 * t48) * t219 + (t18 * t72 + t27 * t51) * t222 + (t27 * t151 - (t51 * t249 + t60 * t72) * t215) * t207 + t72 * t245 + t6 * t234 + (t17 * t68 + t26 * t47) * t220 + (t17 * t71 + t26 * t50) * t223 + (t26 * t149 - (t50 * t250 + t59 * t71) * t216) * t209 + t71 * t246 + t5 * t235 + (t16 * t67 + t25 * t46) * t221 + (t16 * t70 + t25 * t49) * t224 + (t25 * t147 - (t49 * t251 + t58 * t70) * t217) * t211 + t70 * t247 + t4 * t236 + (-g(2) + t139) * m(4); (t48 * t57 + t54 * t69) * t219 + (t51 * t57 + t54 * t72) * t222 + (t57 * t151 - (t151 * t249 - t218) * t206) * t207 + t6 * t202 - (-t60 * t12 - t9 * t249 - t95 * t151 * t191 - 0.2e1 * t42 * (t103 * t129 + (qJ(2,1) * t198 + t80 * t145) * t151 + t145 * qJ(2,1) * t197 + t200) - ((t180 * rSges(3,1) - t248) * t151 - t145 * (t180 * rSges(3,2) + t155)) * m(3)) * t206 + (t47 * t56 + t53 * t68) * t220 + (t50 * t56 + t53 * t71) * t223 + (t56 * t149 - (t149 * t250 - t218) * t208) * t209 + t5 * t203 - (-t59 * t11 - t8 * t250 - t95 * t149 * t193 - 0.2e1 * t41 * (t103 * t128 + (qJ(2,2) * t198 + t80 * t143) * t149 + t143 * qJ(2,2) * t197 + t200) - ((t181 * rSges(3,1) - t248) * t149 - t143 * (t181 * rSges(3,2) + t155)) * m(3)) * t208 + (t46 * t55 + t52 * t67) * t221 + (t49 * t55 + t52 * t70) * t224 + (t55 * t147 - (t147 * t251 - t218) * t210) * t211 + t4 * t204 - (-t58 * t10 - t7 * t251 - t95 * t147 * t195 - 0.2e1 * t40 * (t103 * t127 + (qJ(2,3) * t198 + t80 * t141) * t147 + t141 * qJ(2,3) * t197 + t200) - ((t182 * rSges(3,1) - t248) * t147 - t141 * (t182 * rSges(3,2) + t155)) * m(3)) * t210 + (-g(3) + t138) * m(4);];
tauX  = t1;
