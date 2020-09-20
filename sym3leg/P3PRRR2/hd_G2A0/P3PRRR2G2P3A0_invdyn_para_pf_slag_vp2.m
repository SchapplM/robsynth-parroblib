% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR2G2P3A0
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
%   pkin=[a3,a4]';
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
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:33
% EndTime: 2020-03-09 21:21:35
% DurationCPUTime: 1.66s
% Computational Cost: add. (3276->215), mult. (5985->384), div. (2988->5), fcn. (6126->24), ass. (0->185)
t123 = xDP(2);
t105 = legFrame(1,2);
t94 = cos(t105);
t164 = t94 * t123;
t124 = xDP(1);
t91 = sin(t105);
t167 = t91 * t124;
t208 = t164 + t167;
t104 = legFrame(2,2);
t93 = cos(t104);
t165 = t93 * t123;
t90 = sin(t104);
t168 = t90 * t124;
t207 = t165 + t168;
t103 = legFrame(3,2);
t92 = cos(t103);
t166 = t92 * t123;
t89 = sin(t103);
t169 = t89 * t124;
t206 = t166 + t169;
t205 = mrSges(2,2) * g(3);
t204 = mrSges(3,2) * pkin(1);
t203 = mrSges(3,2) * g(3);
t109 = sin(qJ(3,3));
t115 = cos(qJ(3,3));
t202 = pkin(1) * (mrSges(3,1) * t109 + mrSges(3,2) * t115);
t111 = sin(qJ(3,2));
t117 = cos(qJ(3,2));
t201 = pkin(1) * (mrSges(3,1) * t111 + mrSges(3,2) * t117);
t113 = sin(qJ(3,1));
t119 = cos(qJ(3,1));
t200 = pkin(1) * (mrSges(3,1) * t113 + mrSges(3,2) * t119);
t199 = pkin(1) * t115;
t198 = pkin(1) * t117;
t197 = pkin(1) * t119;
t116 = cos(qJ(2,3));
t126 = 0.1e1 / pkin(2);
t110 = sin(qJ(2,3));
t163 = t109 * t110;
t181 = t126 * ((pkin(2) * t115 + pkin(1)) * t116 - pkin(2) * t163);
t137 = t206 * t181;
t122 = xDP(3);
t100 = qJ(2,3) + qJ(3,3);
t86 = sin(t100);
t172 = t86 * t122;
t128 = 0.1e1 / pkin(1);
t97 = 0.1e1 / t109;
t175 = t128 * t97;
t64 = t115 * t116 - t163;
t193 = t206 * t64 * t175;
t160 = t122 * t128;
t152 = t126 * t160;
t76 = pkin(1) * t110 + pkin(2) * t86;
t55 = t76 * t97 * t152;
t16 = t55 + (-t137 - t172) * t175 + t193;
t19 = -t137 * t175 + t55;
t196 = t16 * t19;
t118 = cos(qJ(2,2));
t112 = sin(qJ(2,2));
t162 = t111 * t112;
t180 = t126 * ((pkin(2) * t117 + pkin(1)) * t118 - pkin(2) * t162);
t136 = t207 * t180;
t101 = qJ(2,2) + qJ(3,2);
t87 = sin(t101);
t171 = t87 * t122;
t98 = 0.1e1 / t111;
t174 = t128 * t98;
t65 = t117 * t118 - t162;
t192 = t207 * t65 * t174;
t77 = pkin(1) * t112 + pkin(2) * t87;
t56 = t77 * t98 * t152;
t17 = t56 + (-t136 - t171) * t174 + t192;
t20 = -t136 * t174 + t56;
t195 = t17 * t20;
t120 = cos(qJ(2,1));
t114 = sin(qJ(2,1));
t161 = t113 * t114;
t179 = t126 * ((pkin(2) * t119 + pkin(1)) * t120 - pkin(2) * t161);
t135 = t208 * t179;
t102 = qJ(2,1) + qJ(3,1);
t88 = sin(t102);
t170 = t88 * t122;
t99 = 0.1e1 / t113;
t173 = t128 * t99;
t66 = t119 * t120 - t161;
t191 = t208 * t66 * t173;
t78 = pkin(1) * t114 + pkin(2) * t88;
t57 = t78 * t99 * t152;
t18 = t57 + (-t135 - t170) * t173 + t191;
t21 = -t135 * t173 + t57;
t194 = t18 * t21;
t121 = mrSges(3,1) * g(3);
t70 = g(1) * t89 + g(2) * t92;
t190 = t86 * (mrSges(3,1) * t70 - t203) + (mrSges(3,2) * t70 + t121) * cos(t100);
t71 = g(1) * t90 + g(2) * t93;
t189 = t87 * (mrSges(3,1) * t71 - t203) + (mrSges(3,2) * t71 + t121) * cos(t101);
t72 = g(1) * t91 + g(2) * t94;
t188 = t88 * (mrSges(3,1) * t72 - t203) + (mrSges(3,2) * t72 + t121) * cos(t102);
t107 = xDDP(2);
t187 = t107 * t92;
t186 = t107 * t93;
t185 = t107 * t94;
t108 = xDDP(1);
t184 = t108 * t89;
t183 = t108 * t90;
t182 = t108 * t91;
t178 = t126 * t76;
t177 = t126 * t77;
t176 = t126 * t78;
t159 = 0.2e1 * pkin(1) * pkin(2);
t158 = mrSges(3,1) * t199;
t157 = mrSges(3,1) * t198;
t156 = mrSges(3,1) * t197;
t127 = pkin(1) ^ 2;
t151 = m(3) * t127 + Ifges(2,3) + Ifges(3,3);
t34 = -t86 * t97 * t160 + t193;
t10 = ((-pkin(2) * t16 - t34 * t199) * t34 - pkin(2) * t196) * t175;
t13 = t55 / 0.2e1 + (-t172 + (-t169 / 0.2e1 - t166 / 0.2e1) * t181) * t175 + t193;
t83 = t109 * t204;
t61 = t151 - 0.2e1 * t83 + 0.2e1 * t158;
t67 = Ifges(3,3) - t83 + t158;
t125 = pkin(2) ^ 2;
t8 = ((t115 * t13 * t159 + t125 * t16 + t127 * t34) * t126 * t34 + (pkin(2) + t199) * t196) * t175;
t95 = m(3) * pkin(1) + mrSges(2,1);
t82 = t95 * g(3);
t1 = -t61 * t10 - t67 * t8 - 0.2e1 * t13 * t19 * t202 + (mrSges(2,2) * t70 + t82) * t116 + (t70 * t95 - t205) * t110 + t190;
t4 = t34 ^ 2 * t202 - Ifges(3,3) * t8 - t10 * t67 + t190;
t150 = t64 * t1 - t4 * t181;
t35 = -t87 * t98 * t160 + t192;
t11 = ((-pkin(2) * t17 - t35 * t198) * t35 - pkin(2) * t195) * t174;
t14 = t56 / 0.2e1 + (-t171 + (-t168 / 0.2e1 - t165 / 0.2e1) * t180) * t174 + t192;
t84 = t111 * t204;
t62 = t151 - 0.2e1 * t84 + 0.2e1 * t157;
t68 = Ifges(3,3) - t84 + t157;
t9 = ((t117 * t14 * t159 + t125 * t17 + t127 * t35) * t126 * t35 + (pkin(2) + t198) * t195) * t174;
t2 = -t62 * t11 - t68 * t9 - 0.2e1 * t14 * t20 * t201 + (mrSges(2,2) * t71 + t82) * t118 + (t71 * t95 - t205) * t112 + t189;
t5 = t35 ^ 2 * t201 - Ifges(3,3) * t9 - t11 * t68 + t189;
t149 = -t5 * t180 + t65 * t2;
t36 = -t88 * t99 * t160 + t191;
t12 = ((-pkin(2) * t18 - t36 * t197) * t36 - pkin(2) * t194) * t173;
t15 = t57 / 0.2e1 + (-t170 + (-t167 / 0.2e1 - t164 / 0.2e1) * t179) * t173 + t191;
t85 = t113 * t204;
t63 = t151 - 0.2e1 * t85 + 0.2e1 * t156;
t69 = Ifges(3,3) - t85 + t156;
t7 = ((t119 * t15 * t159 + t125 * t18 + t127 * t36) * t126 * t36 + (pkin(2) + t197) * t194) * t173;
t3 = -t63 * t12 - t69 * t7 - 0.2e1 * t15 * t21 * t200 + (mrSges(2,2) * t72 + t82) * t120 + (t72 * t95 - t205) * t114 + t188;
t6 = t36 ^ 2 * t200 - Ifges(3,3) * t7 - t12 * t69 + t188;
t148 = -t6 * t179 + t66 * t3;
t133 = (-t67 * t181 + t61 * t64) * t175;
t22 = t89 * t133;
t134 = (-Ifges(3,3) * t181 + t64 * t67) * t175;
t28 = t89 * t134;
t147 = -t28 * t181 + t22 * t64;
t131 = (-t68 * t180 + t62 * t65) * t174;
t23 = t90 * t131;
t132 = (-Ifges(3,3) * t180 + t65 * t68) * t174;
t29 = t90 * t132;
t146 = -t29 * t180 + t23 * t65;
t129 = (-t69 * t179 + t63 * t66) * t173;
t24 = t91 * t129;
t130 = (-Ifges(3,3) * t179 + t66 * t69) * t173;
t30 = t91 * t130;
t145 = -t30 * t179 + t24 * t66;
t25 = t92 * t133;
t31 = t92 * t134;
t144 = -t31 * t181 + t25 * t64;
t26 = t93 * t131;
t32 = t93 * t132;
t143 = -t32 * t180 + t26 * t65;
t27 = t94 * t129;
t33 = t94 * t130;
t142 = -t33 * t179 + t27 * t66;
t138 = -t89 * t92 - t90 * t93 - t91 * t94;
t106 = xDDP(3);
t96 = m(1) + m(2) + m(3);
t75 = g(1) * t94 - g(2) * t91;
t74 = g(1) * t93 - g(2) * t90;
t73 = g(1) * t92 - g(2) * t89;
t42 = (Ifges(3,3) * t176 - t69 * t88) * t173;
t41 = (Ifges(3,3) * t177 - t68 * t87) * t174;
t40 = (Ifges(3,3) * t178 - t67 * t86) * t175;
t39 = (t69 * t176 - t63 * t88) * t173;
t38 = (t68 * t177 - t62 * t87) * t174;
t37 = (t67 * t178 - t61 * t86) * t175;
t43 = [(-g(1) + t108) * m(4) + (-t92 * t73 - t93 * t74 - t94 * t75 + (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) * t108 + t138 * t107) * t96 + (((t30 * t176 - t24 * t88) * t106 + t145 * t185 + (t145 * t108 + t148) * t91) * t99 + ((t29 * t177 - t23 * t87) * t106 + t146 * t186 + (t146 * t108 + t149) * t90) * t98 + ((t28 * t178 - t22 * t86) * t106 + t147 * t187 + (t147 * t108 + t150) * t89) * t97) * t128; (-g(2) + t107) * m(4) + (t89 * t73 + t90 * t74 + t91 * t75 + t138 * t108 + (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) * t107) * t96 + (((t33 * t176 - t27 * t88) * t106 + t142 * t182 + (t142 * t107 + t148) * t94) * t99 + ((t32 * t177 - t26 * t87) * t106 + t143 * t183 + (t143 * t107 + t149) * t93) * t98 + ((t31 * t178 - t25 * t86) * t106 + t144 * t184 + (t144 * t107 + t150) * t92) * t97) * t128; (-g(3) + t106) * m(4) + (((t42 * t176 - t39 * t88) * t106 - t88 * t3 + t6 * t176 + (t185 + t182) * (-t42 * t179 + t39 * t66)) * t99 + ((t41 * t177 - t38 * t87) * t106 - t87 * t2 + t5 * t177 + (t186 + t183) * (-t41 * t180 + t38 * t65)) * t98 + ((t40 * t178 - t37 * t86) * t106 - t86 * t1 + t4 * t178 + (t187 + t184) * (-t40 * t181 + t37 * t64)) * t97) * t128;];
tauX  = t43;
