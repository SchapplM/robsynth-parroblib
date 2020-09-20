% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR1G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:26
% EndTime: 2020-03-09 21:16:28
% DurationCPUTime: 2.17s
% Computational Cost: add. (1953->271), mult. (4386->525), div. (3120->10), fcn. (5355->18), ass. (0->214)
t113 = sin(qJ(3,3));
t127 = xDP(2);
t128 = xDP(1);
t107 = legFrame(3,2);
t86 = sin(t107);
t89 = cos(t107);
t154 = t127 * t86 - t128 * t89;
t129 = 0.1e1 / pkin(2);
t119 = cos(qJ(3,3));
t97 = 0.1e1 / t119;
t186 = t129 * t97;
t55 = t154 * t186;
t52 = t55 ^ 2;
t234 = t113 * t52;
t123 = cos(qJ(3,1));
t102 = t123 ^ 2;
t233 = (0.2e1 * t102 - 0.1e1) * Ifges(3,4);
t121 = cos(qJ(3,2));
t99 = t121 ^ 2;
t232 = (0.2e1 * t99 - 0.1e1) * Ifges(3,4);
t96 = t119 ^ 2;
t231 = (0.2e1 * t96 - 0.1e1) * Ifges(3,4);
t230 = 0.2e1 * Ifges(3,4);
t226 = Ifges(3,5) / 0.2e1;
t225 = -Ifges(3,6) / 0.2e1;
t106 = mrSges(2,2) - mrSges(3,3);
t224 = t106 / 0.2e1;
t223 = t52 * t97;
t116 = sin(qJ(2,2));
t122 = cos(qJ(2,2));
t115 = sin(qJ(3,2));
t156 = mrSges(3,1) * t121 - mrSges(3,2) * t115;
t146 = mrSges(2,1) + t156;
t59 = -t116 * t106 + t146 * t122;
t94 = 0.1e1 / t116;
t222 = t59 * t94;
t118 = sin(qJ(2,1));
t124 = cos(qJ(2,1));
t117 = sin(qJ(3,1));
t155 = mrSges(3,1) * t123 - mrSges(3,2) * t117;
t145 = mrSges(2,1) + t155;
t60 = -t118 * t106 + t145 * t124;
t95 = 0.1e1 / t118;
t221 = t60 * t95;
t76 = Ifges(3,5) * t113 + Ifges(3,6) * t119;
t220 = t76 * t97;
t79 = mrSges(3,1) * t113 + mrSges(3,2) * t119;
t219 = t79 * t97;
t92 = m(1) + m(2) + m(3);
t218 = t92 * t94;
t217 = t92 * t95;
t114 = sin(qJ(2,3));
t93 = 0.1e1 / t114;
t216 = t93 * t97;
t215 = Ifges(3,1) + Ifges(2,3);
t100 = 0.1e1 / t121;
t108 = legFrame(2,2);
t87 = sin(t108);
t90 = cos(t108);
t153 = t127 * t87 - t128 * t90;
t185 = t100 * t129;
t56 = t153 * t185;
t53 = t56 ^ 2;
t214 = t100 * t53;
t213 = t100 * t94;
t103 = 0.1e1 / t123;
t109 = legFrame(1,2);
t88 = sin(t109);
t91 = cos(t109);
t152 = t127 * t88 - t128 * t91;
t184 = t103 * t129;
t57 = t152 * t184;
t54 = t57 ^ 2;
t212 = t103 * t54;
t211 = t103 * t95;
t126 = xDP(3);
t120 = cos(qJ(2,3));
t98 = 0.1e1 / t119 ^ 2;
t170 = t113 * t120 * t98;
t43 = (-t126 * t97 - t154 * t170) * t93 * t129;
t210 = t113 * t43;
t209 = t113 * t55;
t208 = t114 * t79;
t101 = 0.1e1 / t121 ^ 2;
t167 = t101 * t115 * t122;
t44 = (-t100 * t126 - t153 * t167) * t94 * t129;
t207 = t115 * t44;
t206 = t115 * t56;
t80 = mrSges(3,1) * t115 + mrSges(3,2) * t121;
t205 = t116 * t80;
t104 = 0.1e1 / t123 ^ 2;
t166 = t104 * t117 * t124;
t45 = (-t103 * t126 - t152 * t166) * t95 * t129;
t204 = t117 * t45;
t203 = t117 * t57;
t81 = mrSges(3,1) * t117 + mrSges(3,2) * t123;
t202 = t118 * t81;
t40 = t43 ^ 2;
t201 = t119 * t40;
t200 = t120 * t43;
t199 = t120 * t93;
t41 = t44 ^ 2;
t198 = t121 * t41;
t197 = t122 * t44;
t196 = t122 * t94;
t42 = t45 ^ 2;
t195 = t123 * t42;
t194 = t124 * t45;
t193 = t124 * t95;
t192 = t129 * t86;
t191 = t129 * t87;
t190 = t129 * t88;
t189 = t129 * t89;
t188 = t129 * t90;
t187 = t129 * t91;
t183 = t114 * t119;
t182 = t116 * t121;
t181 = t118 * t123;
t61 = -t113 * t89 + t86 * t183;
t180 = t61 * t216;
t64 = t113 * t86 + t89 * t183;
t179 = t64 * t216;
t178 = t92 * t216;
t62 = -t115 * t90 + t87 * t182;
t177 = t62 * t213;
t65 = t115 * t87 + t90 * t182;
t176 = t65 * t213;
t63 = -t117 * t91 + t88 * t181;
t175 = t63 * t211;
t66 = t117 * t88 + t91 * t181;
t174 = t66 * t211;
t173 = t97 * t208;
t172 = t115 * t214;
t171 = t117 * t212;
t169 = t129 * t205;
t168 = t129 * t202;
t165 = t93 * t170;
t164 = t94 * t167;
t163 = t95 * t166;
t70 = g(1) * t86 + g(2) * t89;
t162 = g(3) * t120 + t114 * t70;
t71 = g(1) * t87 + g(2) * t90;
t161 = g(3) * t122 + t116 * t71;
t72 = g(1) * t88 + g(2) * t91;
t160 = g(3) * t124 + t118 * t72;
t159 = t129 * t164;
t158 = t129 * t163;
t157 = mrSges(3,1) * t119 - mrSges(3,2) * t113;
t151 = t87 * t159;
t150 = t90 * t159;
t149 = t88 * t158;
t148 = t91 * t158;
t147 = mrSges(2,1) + t157;
t10 = ((-t114 * t209 + t119 * t200) * t97 * t43 + (t120 * t55 - t183 * t210) * t98 * t55) * t93;
t105 = Ifges(3,1) - Ifges(3,2);
t25 = (-t201 - t223) * t93 * pkin(2);
t125 = mrSges(2,1) * g(3);
t58 = -t114 * t106 + t147 * t120;
t67 = t113 * t119 * t230 - t105 * t96 + t215;
t85 = t106 * g(3);
t4 = -t58 * t25 - t67 * t10 + t220 * t234 + 0.2e1 * t55 * ((t105 * t210 + t55 * t226) * t119 + t209 * t225 + t43 * t231) + (-t147 * t70 + t85) * t120 + t114 * (t157 * g(3) + t106 * t70 + t125);
t73 = g(1) * t89 - g(2) * t86;
t144 = t4 * t165 - t97 * (t25 * t208 - t76 * t10 + (mrSges(3,1) * t73 + t162 * mrSges(3,2)) * t119 + (t162 * mrSges(3,1) - mrSges(3,2) * t73 + Ifges(3,3) * t223 - t105 * t201) * t113 - t40 * t231);
t143 = -Ifges(3,3) * t97 + t76 * t165;
t11 = ((-t116 * t206 + t121 * t197) * t100 * t44 + (t122 * t56 - t182 * t207) * t101 * t56) * t94;
t26 = (-t198 - t214) * t94 * pkin(2);
t68 = t115 * t121 * t230 - t105 * t99 + t215;
t77 = Ifges(3,5) * t115 + Ifges(3,6) * t121;
t5 = -t59 * t26 - t68 * t11 + t77 * t172 + 0.2e1 * t56 * ((t105 * t207 + t56 * t226) * t121 + t206 * t225 + t44 * t232) + (-t146 * t71 + t85) * t122 + t116 * (t156 * g(3) + t106 * t71 + t125);
t74 = g(1) * t90 - g(2) * t87;
t142 = -t100 * (t26 * t205 - t77 * t11 + (mrSges(3,1) * t74 + t161 * mrSges(3,2)) * t121 + (t161 * mrSges(3,1) - mrSges(3,2) * t74 + Ifges(3,3) * t214 - t105 * t198) * t115 - t41 * t232) + t5 * t164;
t12 = ((-t118 * t203 + t123 * t194) * t103 * t45 + (t124 * t57 - t181 * t204) * t104 * t57) * t95;
t27 = (-t195 - t212) * t95 * pkin(2);
t69 = t117 * t123 * t230 - t102 * t105 + t215;
t78 = Ifges(3,5) * t117 + Ifges(3,6) * t123;
t6 = -t60 * t27 - t69 * t12 + t78 * t171 + 0.2e1 * t57 * ((t105 * t204 + t57 * t226) * t123 + t203 * t225 + t45 * t233) + (-t145 * t72 + t85) * t124 + t118 * (t155 * g(3) + t106 * t72 + t125);
t75 = g(1) * t91 - g(2) * t88;
t141 = -t103 * (t27 * t202 - t78 * t12 + (mrSges(3,1) * t75 + t160 * mrSges(3,2)) * t123 + (t160 * mrSges(3,1) - mrSges(3,2) * t75 + Ifges(3,3) * t212 - t105 * t195) * t117 - t42 * t233) + t6 * t163;
t137 = t67 * t165 - t220;
t13 = t137 * t189 + t58 * t180;
t140 = t13 * t165 - (t143 * t189 - t61 * t219) * t97;
t16 = -t137 * t192 + t58 * t179;
t139 = t16 * t165 - (-t143 * t192 - t64 * t219) * t97;
t37 = (t120 * t58 - t67 * t186) * t93;
t138 = t37 * t165 - (-t76 * t93 * t186 - t120 * t79) * t97;
t14 = t68 * t150 + (-t77 * t188 + t62 * t222) * t100;
t136 = -t100 * (t77 * t150 + (-Ifges(3,3) * t188 - t62 * t80) * t100) + t14 * t164;
t17 = -t68 * t151 + (t77 * t191 + t65 * t222) * t100;
t135 = -t100 * (-t77 * t151 + (Ifges(3,3) * t191 - t65 * t80) * t100) + t17 * t164;
t38 = (t122 * t59 - t68 * t185) * t94;
t134 = -t100 * (-t77 * t94 * t185 - t122 * t80) + t38 * t164;
t15 = t69 * t148 + (-t78 * t187 + t63 * t221) * t103;
t133 = -t103 * (t78 * t148 + (-Ifges(3,3) * t187 - t63 * t81) * t103) + t15 * t163;
t18 = -t69 * t149 + (t78 * t190 + t66 * t221) * t103;
t132 = -t103 * (-t78 * t149 + (Ifges(3,3) * t190 - t66 * t81) * t103) + t18 * t163;
t39 = (t124 * t60 - t69 * t184) * t95;
t131 = -t103 * (-t78 * t95 * t184 - t124 * t81) + t39 * t163;
t130 = t58 * t165 + t173;
t112 = xDDP(1);
t111 = xDDP(2);
t110 = xDDP(3);
t48 = (t124 * t92 - t60 * t184) * t95;
t47 = (t122 * t92 - t59 * t185) * t94;
t46 = (t120 * t92 - t58 * t186) * t93;
t24 = -t60 * t149 + (-t88 * t168 + t66 * t217) * t103;
t23 = -t59 * t151 + (-t87 * t169 + t65 * t218) * t100;
t22 = -t130 * t192 + t64 * t178;
t21 = t60 * t148 + (t91 * t168 + t63 * t217) * t103;
t20 = t59 * t150 + (t90 * t169 + t62 * t218) * t100;
t19 = t130 * t189 + t61 * t178;
t3 = -t60 * t12 - t171 * t202 + (-mrSges(2,1) * t42 - t155 * (t42 + t54)) * t118 - 0.2e1 * (t45 * t224 + t81 * t57) * t194 + (-t27 - t72) * t92;
t2 = -t59 * t11 - t172 * t205 + (-mrSges(2,1) * t41 - t156 * (t41 + t53)) * t116 - 0.2e1 * (t44 * t224 + t80 * t56) * t197 + (-t26 - t71) * t92;
t1 = -t58 * t10 - t173 * t234 + (-mrSges(2,1) * t40 - t157 * (t40 + t52)) * t114 - 0.2e1 * (t43 * t224 + t79 * t55) * t200 + (-t25 - t70) * t92;
t7 = [t1 * t180 + t2 * t177 + t3 * t175 - g(1) * m(4) + (t21 * t174 + t20 * t176 + t19 * t179) * t111 + (t19 * t199 + t21 * t193 + t20 * t196) * t110 + (t21 * t175 + t20 * t177 + t19 * t180 + m(4)) * t112 + ((-t133 * t88 - t136 * t87 - t140 * t86) * t111 + (-t13 * t216 - t14 * t213 - t15 * t211) * t110 + (t133 * t112 + t141) * t91 + (t136 * t112 + t142) * t90 + (t140 * t112 + t144) * t89) * t129; t1 * t179 + t2 * t176 + t3 * t174 - g(2) * m(4) + (t24 * t175 + t23 * t177 + t22 * t180) * t112 + (t24 * t193 + t23 * t196 + t22 * t199) * t110 + (t24 * t174 + t23 * t176 + t22 * t179 + m(4)) * t111 + ((t132 * t91 + t135 * t90 + t139 * t89) * t112 + (-t16 * t216 - t17 * t213 - t18 * t211) * t110 + (-t132 * t111 - t141) * t88 + (-t135 * t111 - t142) * t87 + (-t139 * t111 - t144) * t86) * t129; t1 * t199 + t2 * t196 + t3 * t193 - g(3) * m(4) + (t48 * t175 + t47 * t177 + t46 * t180) * t112 + (t48 * t174 + t47 * t176 + t46 * t179) * t111 + (t48 * t193 + t47 * t196 + t46 * t199 + m(4)) * t110 + (-t5 * t213 - t6 * t211 - t4 * t216 + (t131 * t91 + t134 * t90 + t138 * t89) * t112 + (-t131 * t88 - t134 * t87 - t138 * t86) * t111 + (-t39 * t211 - t38 * t213 - t37 * t216) * t110) * t129;];
tauX  = t7;
