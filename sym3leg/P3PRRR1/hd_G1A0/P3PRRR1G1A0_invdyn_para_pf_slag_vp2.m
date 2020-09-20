% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:49
% EndTime: 2020-03-09 21:14:50
% DurationCPUTime: 1.25s
% Computational Cost: add. (9543->177), mult. (5422->348), div. (1542->5), fcn. (6384->24), ass. (0->171)
t110 = pkin(7) + qJ(2,3);
t100 = qJ(3,3) + t110;
t85 = sin(t100);
t88 = cos(t100);
t94 = sin(t110);
t97 = cos(t110);
t200 = 0.1e1 / (t85 * t97 - t94 * t88);
t111 = pkin(7) + qJ(2,2);
t101 = qJ(3,2) + t111;
t86 = sin(t101);
t89 = cos(t101);
t95 = sin(t111);
t98 = cos(t111);
t199 = 0.1e1 / (t86 * t98 - t95 * t89);
t112 = pkin(7) + qJ(2,1);
t102 = qJ(3,1) + t112;
t87 = sin(t102);
t90 = cos(t102);
t96 = sin(t112);
t99 = cos(t112);
t198 = 0.1e1 / (t87 * t99 - t96 * t90);
t197 = mrSges(3,1) * pkin(2);
t196 = mrSges(3,2) * pkin(2);
t120 = sin(qJ(3,3));
t123 = cos(qJ(3,3));
t195 = pkin(2) * (t120 * mrSges(3,1) + t123 * mrSges(3,2));
t121 = sin(qJ(3,2));
t124 = cos(qJ(3,2));
t194 = pkin(2) * (t121 * mrSges(3,1) + t124 * mrSges(3,2));
t122 = sin(qJ(3,1));
t125 = cos(qJ(3,1));
t193 = pkin(2) * (t122 * mrSges(3,1) + t125 * mrSges(3,2));
t129 = 0.1e1 / pkin(3);
t131 = 0.1e1 / pkin(2);
t158 = t129 * t131;
t148 = t200 * t158;
t127 = xDP(1);
t114 = legFrame(3,3);
t103 = sin(t114);
t106 = cos(t114);
t65 = -t85 * t103 + t106 * t88;
t46 = -pkin(2) * (t94 * t103 - t106 * t97) + t65 * pkin(3);
t163 = t46 * t127;
t126 = xDP(2);
t64 = t103 * t88 + t106 * t85;
t43 = pkin(2) * (t103 * t97 + t94 * t106) + t64 * pkin(3);
t166 = t43 * t126;
t134 = (t163 + t166) * t148;
t159 = t127 * t131;
t160 = t126 * t131;
t182 = t200 * t65;
t183 = t200 * t64;
t34 = t159 * t182 + t160 * t183;
t16 = -t134 + t34;
t192 = t16 * t134;
t147 = t199 * t158;
t115 = legFrame(2,3);
t104 = sin(t115);
t107 = cos(t115);
t67 = -t86 * t104 + t107 * t89;
t47 = -pkin(2) * (t95 * t104 - t107 * t98) + t67 * pkin(3);
t162 = t47 * t127;
t66 = t104 * t89 + t107 * t86;
t44 = pkin(2) * (t104 * t98 + t95 * t107) + t66 * pkin(3);
t165 = t44 * t126;
t133 = (t162 + t165) * t147;
t178 = t199 * t67;
t179 = t199 * t66;
t35 = t159 * t178 + t160 * t179;
t17 = -t133 + t35;
t191 = t17 * t133;
t146 = t198 * t158;
t116 = legFrame(1,3);
t105 = sin(t116);
t108 = cos(t116);
t69 = -t87 * t105 + t108 * t90;
t48 = -pkin(2) * (t96 * t105 - t108 * t99) + t69 * pkin(3);
t161 = t48 * t127;
t68 = t105 * t90 + t108 * t87;
t45 = pkin(2) * (t105 * t99 + t96 * t108) + t68 * pkin(3);
t164 = t45 * t126;
t132 = (t161 + t164) * t146;
t174 = t198 * t69;
t175 = t198 * t68;
t36 = t159 * t174 + t160 * t175;
t18 = -t132 + t36;
t190 = t18 * t132;
t189 = t43 * t200;
t188 = t44 * t199;
t187 = t45 * t198;
t186 = t46 * t200;
t185 = t47 * t199;
t184 = t48 * t198;
t130 = pkin(2) ^ 2;
t141 = m(3) * t130 + Ifges(2,3) + Ifges(3,3);
t157 = t123 * t197;
t91 = t120 * t196;
t70 = t141 - 0.2e1 * t91 + 0.2e1 * t157;
t181 = t200 * t70;
t79 = Ifges(3,3) - t91 + t157;
t180 = t200 * t79;
t156 = t124 * t197;
t92 = t121 * t196;
t71 = t141 - 0.2e1 * t92 + 0.2e1 * t156;
t177 = t199 * t71;
t80 = Ifges(3,3) - t92 + t156;
t176 = t199 * t80;
t155 = t125 * t197;
t93 = t122 * t196;
t72 = t141 - 0.2e1 * t93 + 0.2e1 * t155;
t173 = t198 * t72;
t81 = Ifges(3,3) - t93 + t155;
t172 = t198 * t81;
t170 = t129 * t200;
t169 = t129 * t199;
t168 = t129 * t198;
t167 = t131 * t200;
t154 = Ifges(3,3) * t170;
t153 = Ifges(3,3) * t169;
t152 = Ifges(3,3) * t168;
t151 = t79 * t170;
t150 = t80 * t169;
t149 = t81 * t168;
t145 = 0.2e1 * pkin(2) * pkin(3);
t73 = -t103 * g(1) + t106 * g(2);
t76 = t106 * g(1) + t103 * g(2);
t144 = -t88 * (mrSges(3,1) * t73 - mrSges(3,2) * t76) + (mrSges(3,1) * t76 + mrSges(3,2) * t73) * t85;
t74 = -t104 * g(1) + t107 * g(2);
t77 = t107 * g(1) + t104 * g(2);
t143 = -t89 * (mrSges(3,1) * t74 - mrSges(3,2) * t77) + (mrSges(3,1) * t77 + mrSges(3,2) * t74) * t86;
t75 = -t105 * g(1) + t108 * g(2);
t78 = t108 * g(1) + t105 * g(2);
t142 = -t90 * (mrSges(3,1) * t75 - mrSges(3,2) * t78) + (mrSges(3,1) * t78 + mrSges(3,2) * t75) * t87;
t140 = t85 * t94 + t88 * t97;
t139 = t86 * t95 + t89 * t98;
t138 = t87 * t96 + t90 * t99;
t137 = t140 * pkin(2);
t136 = t139 * pkin(2);
t135 = t138 * pkin(2);
t128 = pkin(3) ^ 2;
t119 = xDDP(1);
t118 = xDDP(2);
t109 = m(3) * pkin(2) + mrSges(2,1);
t33 = (-t48 * t152 + t69 * t172) * t131;
t32 = (-t45 * t152 + t68 * t172) * t131;
t31 = (-t47 * t153 + t67 * t176) * t131;
t30 = (-t44 * t153 + t66 * t176) * t131;
t29 = (-t46 * t154 + t65 * t180) * t131;
t28 = (-t43 * t154 + t64 * t180) * t131;
t27 = (-t48 * t149 + t69 * t173) * t131;
t26 = (-t45 * t149 + t68 * t173) * t131;
t25 = (-t47 * t150 + t67 * t177) * t131;
t24 = (-t44 * t150 + t66 * t177) * t131;
t23 = (-t46 * t151 + t65 * t181) * t131;
t22 = (-t43 * t151 + t64 * t181) * t131;
t15 = -(t161 / 0.2e1 + t164 / 0.2e1) * t146 + t36;
t14 = -(t162 / 0.2e1 + t165 / 0.2e1) * t147 + t35;
t13 = -(t163 / 0.2e1 + t166 / 0.2e1) * t148 + t34;
t12 = (-pkin(3) * t191 + (pkin(3) * t17 + t35 * t136) * t35) * t199 * t131;
t11 = (-pkin(3) * t190 + (t18 * pkin(3) + t36 * t135) * t36) * t198 * t131;
t10 = (-pkin(3) * t192 + (pkin(3) * t16 + t34 * t137) * t34) * t167;
t9 = (-(-t138 * t15 * t145 - t18 * t128 - t130 * t36) * t36 * t168 - (pkin(3) + t135) * t198 * t190) * t131;
t8 = (-(-t139 * t14 * t145 - t17 * t128 - t130 * t35) * t35 * t169 - (pkin(3) + t136) * t199 * t191) * t131;
t7 = ((-t140 * t13 * t145 - t16 * t128 - t130 * t34) * t129 * t34 + (pkin(3) + t137) * t192) * t167;
t6 = t36 ^ 2 * t193 - Ifges(3,3) * t9 + t81 * t11 + t142;
t5 = t35 ^ 2 * t194 - Ifges(3,3) * t8 + t80 * t12 + t143;
t4 = t34 ^ 2 * t195 + Ifges(3,3) * t7 + t79 * t10 + t144;
t3 = t72 * t11 - t81 * t9 + 0.2e1 * t15 * t132 * t193 + (mrSges(2,2) * t78 - t75 * t109) * t99 + t96 * (t75 * mrSges(2,2) + t109 * t78) + t142;
t2 = t71 * t12 - t80 * t8 + 0.2e1 * t14 * t133 * t194 + (mrSges(2,2) * t77 - t74 * t109) * t98 + t95 * (t74 * mrSges(2,2) + t109 * t77) + t143;
t1 = t70 * t10 + t79 * t7 + 0.2e1 * t13 * t134 * t195 + (mrSges(2,2) * t76 - t73 * t109) * t97 + t94 * (t73 * mrSges(2,2) + t109 * t76) + t144;
t19 = [(-g(1) + t119) * m(4) + (t1 * t182 + t2 * t178 + t3 * t174 + (-t6 * t184 - t5 * t185 - t4 * t186) * t129 + (t23 * t182 + t25 * t178 + t27 * t174 + (-t33 * t184 - t31 * t185 - t29 * t186) * t129) * t119 + (t23 * t183 + t25 * t179 + t27 * t175 + (-t33 * t187 - t31 * t188 - t29 * t189) * t129) * t118) * t131; (-g(2) + t118) * m(4) + (t1 * t183 + t2 * t179 + t3 * t175 + (-t6 * t187 - t5 * t188 - t4 * t189) * t129 + (t22 * t182 + t24 * t178 + t26 * t174 + (-t32 * t184 - t30 * t185 - t28 * t186) * t129) * t119 + (t22 * t183 + t24 * t179 + t26 * t175 + (-t32 * t187 - t30 * t188 - t28 * t189) * t129) * t118) * t131; (-g(3) + xDDP(3)) * ((3 * m(1)) + (3 * m(2)) + 0.3e1 * m(3) + m(4));];
tauX  = t19;
