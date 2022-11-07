% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR1V1G2A0
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
%   pkin=[a3,a4,d4]';
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
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:07:56
% EndTime: 2022-11-04 17:07:58
% DurationCPUTime: 2.12s
% Computational Cost: add. (8169->347), mult. (9480->635), div. (2637->7), fcn. (7038->24), ass. (0->225)
t137 = legFrame(3,2);
t109 = sin(t137);
t112 = cos(t137);
t162 = xDP(2);
t163 = xDP(1);
t149 = cos(qJ(2,3));
t124 = 0.1e1 / t149;
t164 = pkin(2) + pkin(1);
t130 = 0.1e1 / t164;
t216 = t124 * t130;
t67 = (t109 * t163 + t112 * t162) * t216;
t64 = t67 ^ 2;
t138 = legFrame(2,2);
t110 = sin(t138);
t113 = cos(t138);
t151 = cos(qJ(2,2));
t126 = 0.1e1 / t151;
t215 = t126 * t130;
t68 = (t110 * t163 + t113 * t162) * t215;
t65 = t68 ^ 2;
t139 = legFrame(1,2);
t111 = sin(t139);
t114 = cos(t139);
t153 = cos(qJ(2,1));
t128 = 0.1e1 / t153;
t214 = t128 * t130;
t69 = (t111 * t163 + t114 * t162) * t214;
t66 = t69 ^ 2;
t134 = pkin(3) + qJ(3,3);
t119 = 0.1e1 / t134;
t150 = cos(qJ(1,3));
t161 = xDP(3);
t143 = sin(qJ(2,3));
t144 = sin(qJ(1,3));
t212 = t144 * t149;
t76 = -t109 * t212 + t143 * t112;
t79 = t109 * t143 + t112 * t212;
t43 = (t150 * t161 + (t162 * t76 + t163 * t79) * t124) * t119;
t256 = 0.2e1 * t43;
t135 = pkin(3) + qJ(3,2);
t120 = 0.1e1 / t135;
t152 = cos(qJ(1,2));
t145 = sin(qJ(2,2));
t146 = sin(qJ(1,2));
t210 = t146 * t151;
t77 = -t110 * t210 + t145 * t113;
t80 = t110 * t145 + t113 * t210;
t44 = (t152 * t161 + (t162 * t77 + t163 * t80) * t126) * t120;
t255 = 0.2e1 * t44;
t136 = pkin(3) + qJ(3,1);
t121 = 0.1e1 / t136;
t154 = cos(qJ(1,1));
t147 = sin(qJ(2,1));
t148 = sin(qJ(1,1));
t208 = t148 * t153;
t78 = -t111 * t208 + t147 * t114;
t81 = t111 * t147 + t114 * t208;
t45 = (t154 * t161 + (t162 * t78 + t163 * t81) * t128) * t121;
t254 = 0.2e1 * t45;
t129 = t164 ^ 2;
t253 = m(3) / 0.2e1;
t155 = pkin(1) + rSges(3,1);
t169 = rSges(2,2) ^ 2;
t171 = rSges(2,1) ^ 2;
t174 = (rSges(3,2) + t155) * (-rSges(3,2) + t155) * m(3) - Icges(2,1) - Icges(3,1) + Icges(2,2) + Icges(3,2) + (-t169 + t171) * m(2);
t252 = t174 / 0.2e1;
t251 = m(1) * rSges(1,1);
t250 = m(2) * rSges(2,2);
t98 = -t143 * rSges(3,2) + t155 * t149;
t249 = m(3) * t98;
t99 = -t145 * rSges(3,2) + t155 * t151;
t248 = m(3) * t99;
t247 = rSges(2,3) * m(2);
t100 = -t147 * rSges(3,2) + t155 * t153;
t246 = m(3) * t100;
t245 = m(3) * t119;
t244 = m(3) * t120;
t243 = m(3) * t121;
t131 = rSges(3,3) + qJ(3,3);
t242 = t131 * m(3);
t132 = rSges(3,3) + qJ(3,2);
t241 = t132 * m(3);
t133 = rSges(3,3) + qJ(3,1);
t240 = t133 * m(3);
t239 = t155 * m(3);
t97 = rSges(2,1) * t250 + rSges(3,2) * t239 - Icges(2,4) - Icges(3,4);
t238 = t97 * t67;
t237 = t97 * t68;
t236 = t97 * t69;
t193 = rSges(2,2) * t247 - Icges(2,6) - Icges(3,6);
t177 = rSges(3,2) * t242 + t193;
t181 = rSges(2,1) * t247 - Icges(2,5) - Icges(3,5);
t88 = t131 * t239 + t181;
t55 = -t88 * t143 - t149 * t177;
t235 = t119 * t55;
t176 = rSges(3,2) * t241 + t193;
t89 = t132 * t239 + t181;
t56 = -t89 * t145 - t151 * t176;
t234 = t120 * t56;
t175 = rSges(3,2) * t240 + t193;
t90 = t133 * t239 + t181;
t57 = -t90 * t147 - t153 * t175;
t233 = t121 * t57;
t232 = t124 * t76;
t231 = t124 * t79;
t230 = t124 * t98;
t229 = t126 * t77;
t228 = t126 * t80;
t227 = t126 * t99;
t226 = t128 * t78;
t225 = t128 * t81;
t191 = rSges(3,2) ^ 2 + (pkin(1) ^ 2) + ((2 * pkin(1) + rSges(3,1)) * rSges(3,1));
t201 = t169 + t171;
t75 = t201 * m(2) + t191 * m(3) + Icges(2,3) + Icges(3,3);
t224 = t130 * t75;
t223 = t149 * t174;
t222 = t151 * t174;
t221 = t153 * t174;
t220 = t64 * t124;
t219 = t65 * t126;
t218 = t66 * t128;
t217 = t100 * t128;
t213 = t143 * t164;
t211 = t145 * t164;
t209 = t147 * t164;
t106 = m(2) * rSges(2,1) + t239;
t207 = t149 * t106;
t206 = t151 * t106;
t205 = t153 * t106;
t204 = t164 * t149;
t203 = t164 * t151;
t202 = t164 * t153;
t200 = 0.2e1 * t164;
t180 = -t134 * t150 + t144 * t204;
t58 = -t180 * t109 + t112 * t213;
t61 = t109 * t213 + t180 * t112;
t85 = t144 * t134 + t150 * t204;
t25 = (t161 * t85 + t162 * t58 + t163 * t61) * t119;
t179 = -t135 * t152 + t146 * t203;
t59 = -t179 * t110 + t113 * t211;
t62 = t110 * t211 + t179 * t113;
t86 = t146 * t135 + t152 * t203;
t26 = (t161 * t86 + t162 * t59 + t163 * t62) * t120;
t178 = -t136 * t154 + t148 * t202;
t60 = -t178 * t111 + t114 * t209;
t63 = t111 * t209 + t178 * t114;
t87 = t148 * t136 + t154 * t202;
t27 = (t161 * t87 + t162 * t60 + t163 * t63) * t121;
t199 = t55 * t216;
t198 = t56 * t215;
t197 = t57 * t214;
t196 = t130 * t150 * t55;
t195 = t130 * t152 * t56;
t194 = t130 * t154 * t57;
t192 = m(1) * rSges(1,2) - t247;
t94 = t112 * g(1) - t109 * g(2);
t190 = g(3) * t150 + t144 * t94;
t189 = g(3) * t144 - t94 * t150;
t95 = t113 * g(1) - t110 * g(2);
t188 = g(3) * t152 + t146 * t95;
t187 = g(3) * t146 - t95 * t152;
t96 = t114 * g(1) - t111 * g(2);
t186 = g(3) * t154 + t148 * t96;
t185 = g(3) * t148 - t96 * t154;
t184 = rSges(3,2) * t149 + t155 * t143;
t183 = rSges(3,2) * t151 + t155 * t145;
t182 = rSges(3,2) * t153 + t155 * t147;
t173 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(2,3) ^ 2 + t201) * m(2) / 0.2e1;
t167 = 0.2e1 * qJ(2,1);
t166 = 0.2e1 * qJ(2,2);
t165 = 0.2e1 * qJ(2,3);
t142 = xDDP(1);
t141 = xDDP(2);
t140 = xDDP(3);
t127 = t153 ^ 2;
t125 = t151 ^ 2;
t123 = t149 ^ 2;
t118 = g(3) * t251;
t107 = m(3) * rSges(3,2) + t250;
t103 = t192 - t240;
t102 = t192 - t241;
t101 = t192 - t242;
t93 = t111 * g(1) + t114 * g(2);
t92 = t110 * g(1) + t113 * g(2);
t91 = t109 * g(1) + t112 * g(2);
t42 = (-t81 * t217 + t63) * t243;
t41 = (-t80 * t227 + t62) * t244;
t40 = (-t79 * t230 + t61) * t245;
t39 = (-t78 * t217 + t60) * t243;
t38 = (-t77 * t227 + t59) * t244;
t37 = (-t76 * t230 + t58) * t245;
t36 = cos(t167) * t252 + (0.2e1 * t133 ^ 2 + t191) * t253 - t97 * sin(t167) + t173;
t35 = cos(t166) * t252 + (0.2e1 * t132 ^ 2 + t191) * t253 - t97 * sin(t166) + t173;
t34 = cos(t165) * t252 + (0.2e1 * t131 ^ 2 + t191) * t253 - t97 * sin(t165) + t173;
t33 = (t111 * t224 + t81 * t233) * t128;
t32 = (t110 * t224 + t80 * t234) * t126;
t31 = (t109 * t224 + t79 * t235) * t124;
t30 = (t114 * t224 + t78 * t233) * t128;
t29 = (t113 * t224 + t77 * t234) * t126;
t28 = (t112 * t224 + t76 * t235) * t124;
t24 = (t154 * t36 - t87 * t246) * t121;
t23 = (t152 * t35 - t86 * t248) * t120;
t22 = (t150 * t34 - t85 * t249) * t119;
t21 = t111 * t197 + (t36 * t225 - t63 * t246) * t121;
t20 = t110 * t198 + (t35 * t228 - t62 * t248) * t120;
t19 = t109 * t199 + (t34 * t231 - t61 * t249) * t119;
t18 = t114 * t197 + (t36 * t226 - t60 * t246) * t121;
t17 = t113 * t198 + (t35 * t229 - t59 * t248) * t120;
t16 = t112 * t199 + (t34 * t232 - t58 * t249) * t119;
t15 = (-t164 * t218 + (-t45 * t202 + 0.2e1 * t27) * t45) * t121;
t14 = (-t164 * t219 + (-t44 * t203 + 0.2e1 * t26) * t44) * t120;
t13 = (-t164 * t220 + (-t43 * t204 + 0.2e1 * t25) * t43) * t119;
t12 = (-t129 * t66 + ((-t129 * t127 - t136 ^ 2) * t45 + (t147 * t69 * t136 + t27 * t153) * t200) * t45) * t121;
t11 = (-t129 * t65 + ((-t129 * t125 - t135 ^ 2) * t44 + (t145 * t68 * t135 + t26 * t151) * t200) * t44) * t120;
t10 = (-t129 * t64 + ((-t129 * t123 - t134 ^ 2) * t43 + (t143 * t67 * t134 + t25 * t149) * t200) * t43) * t119;
t9 = -t57 * t15 + t75 * t147 * t218 + (-t182 * t27 * m(3) + (t147 * t221 / 0.2e1 + (t127 - 0.1e1 / 0.2e1) * t97) * t45) * t254 + (t186 * t106 + t93 * t107) * t147 + t153 * (-t93 * t106 + t186 * t107);
t8 = -t56 * t14 + t75 * t145 * t219 + (-t183 * t26 * m(3) + (t145 * t222 / 0.2e1 + (t125 - 0.1e1 / 0.2e1) * t97) * t44) * t255 + (t188 * t106 + t92 * t107) * t145 + t151 * (-t92 * t106 + t188 * t107);
t7 = -t55 * t13 + t75 * t143 * t220 + (-t184 * t25 * m(3) + (t143 * t223 / 0.2e1 + (t123 - 0.1e1 / 0.2e1) * t97) * t43) * t256 + (t190 * t106 + t91 * t107) * t143 + t149 * (-t91 * t106 + t190 * t107);
t6 = t15 * t246 + (-t12 + (-t45 * t133 / 0.2e1 + t182 * t69) * t254 - t185) * m(3);
t5 = t14 * t248 + (-t11 + (-t44 * t132 / 0.2e1 + t183 * t68) * t255 - t187) * m(3);
t4 = t13 * t249 + (-t10 + (-t43 * t131 / 0.2e1 + t184 * t67) * t256 - t189) * m(3);
t3 = -t36 * t15 + t12 * t246 - 0.4e1 * t45 * t127 * t236 - t66 * t90 * t153 + (t27 * t240 + t236) * t254 + (t103 * g(3) + (-t205 - t251) * t96) * t154 + (g(3) * t205 + t96 * t103 + t118) * t148 + (-0.2e1 * t69 * t45 * t221 + (t57 * t128 + t175) * t66 - t185 * t107) * t147;
t2 = -t35 * t14 + t11 * t248 - 0.4e1 * t44 * t125 * t237 - t65 * t89 * t151 + (t26 * t241 + t237) * t255 + (t102 * g(3) + (-t206 - t251) * t95) * t152 + (g(3) * t206 + t95 * t102 + t118) * t146 + (-0.2e1 * t68 * t44 * t222 + (t56 * t126 + t176) * t65 - t187 * t107) * t145;
t1 = -t34 * t13 + t10 * t249 - 0.4e1 * t43 * t123 * t238 - t64 * t88 * t149 + (t25 * t242 + t238) * t256 + (t101 * g(3) + (-t207 - t251) * t94) * t150 + (g(3) * t207 + t94 * t101 + t118) * t144 + (-0.2e1 * t67 * t43 * t223 + (t55 * t124 + t177) * t64 - t189 * t107) * t143;
t46 = [(-g(1) + t142) * m(4) + ((t21 * t225 + t42 * t63) * t142 + (t21 * t226 + t42 * t60) * t141 + (t154 * t21 + t42 * t87) * t140 + t3 * t225 + t63 * t6) * t121 + ((t114 * t141 * t33 + (t142 * t33 + t9) * t111) * t128 + (t113 * t141 * t32 + (t142 * t32 + t8) * t110) * t126 + (t112 * t141 * t31 + (t142 * t31 + t7) * t109) * t124) * t130 + ((t20 * t228 + t41 * t62) * t142 + (t20 * t229 + t41 * t59) * t141 + (t152 * t20 + t41 * t86) * t140 + t2 * t228 + t62 * t5) * t120 + ((t19 * t231 + t40 * t61) * t142 + (t19 * t232 + t40 * t58) * t141 + (t150 * t19 + t40 * t85) * t140 + t1 * t231 + t61 * t4) * t119; (-g(2) + t141) * m(4) + ((t18 * t225 + t39 * t63) * t142 + (t18 * t226 + t39 * t60) * t141 + (t154 * t18 + t39 * t87) * t140 + t3 * t226 + t60 * t6) * t121 + ((t111 * t142 * t30 + (t141 * t30 + t9) * t114) * t128 + (t110 * t142 * t29 + (t141 * t29 + t8) * t113) * t126 + (t109 * t142 * t28 + (t141 * t28 + t7) * t112) * t124) * t130 + ((t17 * t228 + t38 * t62) * t142 + (t17 * t229 + t38 * t59) * t141 + (t152 * t17 + t38 * t86) * t140 + t2 * t229 + t59 * t5) * t120 + ((t16 * t231 + t37 * t61) * t142 + (t16 * t232 + t37 * t58) * t141 + (t150 * t16 + t37 * t85) * t140 + t1 * t232 + t58 * t4) * t119; (-g(3) + t140) * m(4) + (t87 * t6 + (t24 * t140 + t3) * t154 + (t87 * t140 + t60 * t141 + t63 * t142) * (-t100 * t154 + t87) * t243 + ((t111 * t194 + t24 * t81) * t142 + (t114 * t194 + t24 * t78) * t141) * t128) * t121 + (t86 * t5 + (t23 * t140 + t2) * t152 + (t86 * t140 + t59 * t141 + t62 * t142) * (-t152 * t99 + t86) * t244 + ((t110 * t195 + t23 * t80) * t142 + (t113 * t195 + t23 * t77) * t141) * t126) * t120 + (t85 * t4 + (t22 * t140 + t1) * t150 + (t85 * t140 + t58 * t141 + t61 * t142) * (-t150 * t98 + t85) * t245 + ((t109 * t196 + t22 * t79) * t142 + (t112 * t196 + t22 * t76) * t141) * t124) * t119;];
tauX  = t46;
