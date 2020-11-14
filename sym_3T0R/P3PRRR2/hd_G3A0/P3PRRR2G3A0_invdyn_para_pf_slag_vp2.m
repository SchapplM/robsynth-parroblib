% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR2G3A0
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
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR2G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:03
% EndTime: 2020-03-09 21:20:04
% DurationCPUTime: 0.84s
% Computational Cost: add. (2283->162), mult. (1900->295), div. (1374->5), fcn. (1662->36), ass. (0->141)
t155 = mrSges(3,2) * pkin(1);
t108 = sin(qJ(3,3));
t111 = cos(qJ(3,3));
t154 = pkin(1) * (mrSges(3,1) * t108 + mrSges(3,2) * t111);
t109 = sin(qJ(3,2));
t112 = cos(qJ(3,2));
t153 = pkin(1) * (mrSges(3,1) * t109 + mrSges(3,2) * t112);
t110 = sin(qJ(3,1));
t113 = cos(qJ(3,1));
t152 = pkin(1) * (mrSges(3,1) * t110 + mrSges(3,2) * t113);
t150 = pkin(1) * t111;
t149 = pkin(1) * t112;
t148 = pkin(1) * t113;
t115 = xDP(1);
t117 = 0.1e1 / pkin(2);
t119 = 0.1e1 / pkin(1);
t127 = t117 * t119;
t121 = t115 * t127;
t102 = legFrame(3,2);
t85 = -t102 + qJ(2,3);
t79 = qJ(3,3) + t85;
t73 = sin(t79);
t49 = pkin(2) * t73 + pkin(1) * sin(t85);
t96 = 0.1e1 / t108;
t43 = t49 * t96 * t121;
t114 = xDP(2);
t122 = t114 * t127;
t76 = cos(t79);
t52 = -pkin(2) * t76 - pkin(1) * cos(t85);
t46 = t52 * t96 * t122;
t25 = t43 + t46;
t128 = t115 * t119;
t129 = t114 * t119;
t40 = (-t128 * t73 + t129 * t76) * t96;
t16 = t40 + t25;
t147 = t16 * t25;
t103 = legFrame(2,2);
t86 = -t103 + qJ(2,2);
t80 = qJ(3,2) + t86;
t74 = sin(t80);
t50 = pkin(2) * t74 + pkin(1) * sin(t86);
t97 = 0.1e1 / t109;
t44 = t50 * t97 * t121;
t77 = cos(t80);
t53 = -pkin(2) * t77 - pkin(1) * cos(t86);
t47 = t53 * t97 * t122;
t26 = t44 + t47;
t41 = (-t128 * t74 + t129 * t77) * t97;
t17 = t41 + t26;
t146 = t17 * t26;
t104 = legFrame(1,2);
t87 = -t104 + qJ(2,1);
t81 = qJ(3,1) + t87;
t75 = sin(t81);
t51 = pkin(2) * t75 + pkin(1) * sin(t87);
t98 = 0.1e1 / t110;
t45 = t51 * t98 * t121;
t78 = cos(t81);
t54 = -pkin(2) * t78 - pkin(1) * cos(t87);
t48 = t54 * t98 * t122;
t27 = t45 + t48;
t42 = (-t128 * t75 + t129 * t78) * t98;
t18 = t42 + t27;
t145 = t18 * t27;
t88 = sin(t102);
t91 = cos(t102);
t64 = g(1) * t88 + g(2) * t91;
t67 = g(1) * t91 - g(2) * t88;
t99 = qJ(2,3) + qJ(3,3);
t144 = sin(t99) * (mrSges(3,1) * t64 - mrSges(3,2) * t67) + (mrSges(3,1) * t67 + mrSges(3,2) * t64) * cos(t99);
t100 = qJ(2,2) + qJ(3,2);
t89 = sin(t103);
t92 = cos(t103);
t65 = g(1) * t89 + g(2) * t92;
t68 = g(1) * t92 - g(2) * t89;
t143 = sin(t100) * (mrSges(3,1) * t65 - mrSges(3,2) * t68) + (mrSges(3,1) * t68 + mrSges(3,2) * t65) * cos(t100);
t101 = qJ(2,1) + qJ(3,1);
t90 = sin(t104);
t93 = cos(t104);
t66 = g(1) * t90 + g(2) * t93;
t69 = g(1) * t93 - g(2) * t90;
t142 = sin(t101) * (mrSges(3,1) * t66 - mrSges(3,2) * t69) + (mrSges(3,1) * t69 + mrSges(3,2) * t66) * cos(t101);
t141 = t117 * t49;
t140 = t117 * t50;
t139 = t117 * t51;
t138 = t117 * t52;
t137 = t117 * t53;
t136 = t117 * t54;
t125 = mrSges(3,1) * t150;
t82 = t108 * t155;
t61 = Ifges(3,3) - t82 + t125;
t135 = t117 * t61;
t124 = mrSges(3,1) * t149;
t83 = t109 * t155;
t62 = Ifges(3,3) - t83 + t124;
t134 = t117 * t62;
t123 = mrSges(3,1) * t148;
t84 = t110 * t155;
t63 = Ifges(3,3) - t84 + t123;
t133 = t117 * t63;
t132 = t119 * t96;
t131 = t119 * t97;
t130 = t119 * t98;
t126 = 0.2e1 * pkin(1) * pkin(2);
t118 = pkin(1) ^ 2;
t120 = m(3) * t118 + Ifges(2,3) + Ifges(3,3);
t116 = pkin(2) ^ 2;
t107 = xDDP(1);
t106 = xDDP(2);
t94 = m(3) * pkin(1) + mrSges(2,1);
t57 = t120 - 0.2e1 * t84 + 0.2e1 * t123;
t56 = t120 - 0.2e1 * t83 + 0.2e1 * t124;
t55 = t120 - 0.2e1 * t82 + 0.2e1 * t125;
t33 = (Ifges(3,3) * t136 + t63 * t78) * t130;
t32 = (Ifges(3,3) * t137 + t62 * t77) * t131;
t31 = (Ifges(3,3) * t138 + t61 * t76) * t132;
t30 = (Ifges(3,3) * t139 - t63 * t75) * t130;
t29 = (Ifges(3,3) * t140 - t62 * t74) * t131;
t28 = (Ifges(3,3) * t141 - t61 * t73) * t132;
t24 = (t54 * t133 + t57 * t78) * t130;
t23 = (t53 * t134 + t56 * t77) * t131;
t22 = (t52 * t135 + t55 * t76) * t132;
t21 = (t51 * t133 - t57 * t75) * t130;
t20 = (t50 * t134 - t56 * t74) * t131;
t19 = (t135 * t49 - t55 * t73) * t132;
t15 = t45 / 0.2e1 + t48 / 0.2e1 + t42;
t14 = t44 / 0.2e1 + t47 / 0.2e1 + t41;
t13 = t43 / 0.2e1 + t46 / 0.2e1 + t40;
t12 = ((-pkin(2) * t18 - t42 * t148) * t42 - pkin(2) * t145) * t130;
t11 = ((-pkin(2) * t17 - t41 * t149) * t41 - pkin(2) * t146) * t131;
t10 = ((-pkin(2) * t16 - t40 * t150) * t40 - pkin(2) * t147) * t132;
t9 = ((t112 * t126 * t14 + t116 * t17 + t118 * t41) * t117 * t41 + (pkin(2) + t149) * t146) * t131;
t8 = ((t111 * t126 * t13 + t116 * t16 + t118 * t40) * t117 * t40 + (pkin(2) + t150) * t147) * t132;
t7 = ((t113 * t126 * t15 + t116 * t18 + t118 * t42) * t117 * t42 + (pkin(2) + t148) * t145) * t130;
t6 = t42 ^ 2 * t152 - Ifges(3,3) * t7 - t12 * t63 + t142;
t5 = t41 ^ 2 * t153 - Ifges(3,3) * t9 - t11 * t62 + t143;
t4 = t40 ^ 2 * t154 - Ifges(3,3) * t8 - t10 * t61 + t144;
t3 = -t57 * t12 - t63 * t7 - 0.2e1 * t27 * t15 * t152 + (mrSges(2,2) * t66 + t69 * t94) * cos(qJ(2,1)) + sin(qJ(2,1)) * (-mrSges(2,2) * t69 + t66 * t94) + t142;
t2 = -t56 * t11 - t62 * t9 - 0.2e1 * t26 * t14 * t153 + (mrSges(2,2) * t65 + t68 * t94) * cos(qJ(2,2)) + sin(qJ(2,2)) * (-mrSges(2,2) * t68 + t65 * t94) + t143;
t1 = -t55 * t10 - t61 * t8 - 0.2e1 * t25 * t13 * t154 + (mrSges(2,2) * t64 + t67 * t94) * cos(qJ(2,3)) + sin(qJ(2,3)) * (-mrSges(2,2) * t67 + t64 * t94) + t144;
t34 = [(-g(1) + t107) * m(4) + (((t30 * t139 - t21 * t75) * t107 + (t136 * t30 + t21 * t78) * t106 - t75 * t3 + t6 * t139) * t98 + ((t29 * t140 - t20 * t74) * t107 + (t137 * t29 + t20 * t77) * t106 - t74 * t2 + t5 * t140) * t97 + ((t28 * t141 - t19 * t73) * t107 + (t138 * t28 + t19 * t76) * t106 - t73 * t1 + t4 * t141) * t96) * t119; (-g(2) + t106) * m(4) + (((t33 * t139 - t24 * t75) * t107 + (t136 * t33 + t24 * t78) * t106 + t78 * t3 + t6 * t136) * t98 + ((t32 * t140 - t23 * t74) * t107 + (t137 * t32 + t23 * t77) * t106 + t77 * t2 + t5 * t137) * t97 + ((t31 * t141 - t22 * t73) * t107 + (t138 * t31 + t22 * t76) * t106 + t76 * t1 + t4 * t138) * t96) * t119; (-g(3) + xDDP(3)) * ((3 * m(1)) + (3 * m(2)) + 0.3e1 * m(3) + m(4));];
tauX  = t34;
