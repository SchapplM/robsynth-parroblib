% Calculate vector of inverse dynamics forces for parallel robot
% P3PRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRR1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:36
% EndTime: 2019-05-03 14:47:37
% DurationCPUTime: 0.62s
% Computational Cost: add. (789->162), mult. (1402->297), div. (450->4), fcn. (1400->14), ass. (0->114)
t96 = m(1) + m(2);
t116 = pkin(2) * t96;
t106 = 0.1e1 / pkin(2);
t81 = legFrame(3,3);
t66 = sin(t81);
t115 = t106 * t66;
t82 = legFrame(2,3);
t67 = sin(t82);
t114 = t106 * t67;
t83 = legFrame(1,3);
t68 = sin(t83);
t113 = t106 * t68;
t69 = cos(t81);
t112 = t106 * t69;
t70 = cos(t82);
t111 = t106 * t70;
t71 = cos(t83);
t110 = t106 * t71;
t87 = sin(qJ(2,3));
t77 = 0.1e1 / t87;
t109 = t106 * t77;
t88 = sin(qJ(2,2));
t78 = 0.1e1 / t88;
t108 = t106 * t78;
t89 = sin(qJ(2,1));
t79 = 0.1e1 / t89;
t107 = t106 * t79;
t105 = koppelP(1,1);
t104 = koppelP(2,1);
t103 = koppelP(3,1);
t102 = koppelP(1,2);
t101 = koppelP(2,2);
t100 = koppelP(3,2);
t99 = mrSges(3,1);
t98 = mrSges(3,2);
t97 = xP(3);
t95 = xDP(1);
t94 = xDP(2);
t93 = xDP(3);
t92 = cos(qJ(2,1));
t91 = cos(qJ(2,2));
t90 = cos(qJ(2,3));
t86 = xDDP(1);
t85 = xDDP(2);
t84 = xDDP(3);
t80 = t93 ^ 2;
t76 = cos(t97);
t75 = sin(t97);
t65 = mrSges(2,1) * t92 - t89 * mrSges(2,2);
t64 = mrSges(2,1) * t91 - t88 * mrSges(2,2);
t63 = mrSges(2,1) * t90 - t87 * mrSges(2,2);
t62 = g(1) * t71 + g(2) * t68;
t61 = g(1) * t70 + g(2) * t67;
t60 = g(1) * t69 + g(2) * t66;
t59 = -g(1) * t68 + g(2) * t71;
t58 = -g(1) * t67 + g(2) * t70;
t57 = -g(1) * t66 + g(2) * t69;
t56 = -t102 * t75 + t105 * t76;
t55 = -t101 * t75 + t104 * t76;
t54 = -t100 * t75 + t103 * t76;
t53 = t102 * t76 + t105 * t75;
t52 = t101 * t76 + t104 * t75;
t51 = t100 * t76 + t103 * t75;
t50 = -t75 * t98 + t76 * t99;
t49 = t75 * t99 + t76 * t98;
t48 = t68 * t92 + t71 * t89;
t47 = -t68 * t89 + t71 * t92;
t46 = t67 * t91 + t70 * t88;
t45 = -t67 * t88 + t70 * t91;
t44 = t66 * t90 + t69 * t87;
t43 = -t66 * t87 + t69 * t90;
t42 = -t53 * t84 - t56 * t80 + t86;
t41 = -t52 * t84 - t55 * t80 + t86;
t40 = -t51 * t84 - t54 * t80 + t86;
t39 = -t53 * t80 + t56 * t84 + t85;
t38 = -t52 * t80 + t55 * t84 + t85;
t37 = -t51 * t80 + t54 * t84 + t85;
t36 = (-Ifges(2,3) * t113 + t48 * t65) * t79;
t35 = (-Ifges(2,3) * t110 + t47 * t65) * t79;
t34 = (-Ifges(2,3) * t114 + t46 * t64) * t78;
t33 = (-Ifges(2,3) * t111 + t45 * t64) * t78;
t32 = (-Ifges(2,3) * t115 + t44 * t63) * t77;
t31 = (-Ifges(2,3) * t112 + t43 * t63) * t77;
t30 = (-t65 * t113 + t48 * t96) * t79;
t29 = (-t65 * t110 + t47 * t96) * t79;
t28 = (-t64 * t114 + t46 * t96) * t78;
t27 = (-t64 * t111 + t45 * t96) * t78;
t26 = (-t63 * t115 + t44 * t96) * t77;
t25 = (-t63 * t112 + t43 * t96) * t77;
t24 = (t53 * t71 - t56 * t68) * t107;
t23 = (t52 * t70 - t55 * t67) * t108;
t22 = (t51 * t69 - t54 * t66) * t109;
t21 = (-t71 * (-t53 * t93 + t95) - t68 * (t56 * t93 + t94)) * t107;
t20 = (-t70 * (-t52 * t93 + t95) - t67 * (t55 * t93 + t94)) * t108;
t19 = (-t69 * (-t51 * t93 + t95) - t66 * (t54 * t93 + t94)) * t109;
t18 = t21 ^ 2;
t17 = t20 ^ 2;
t16 = t19 ^ 2;
t15 = (-t47 * t53 + t48 * t56) * t79;
t14 = (-t45 * t52 + t46 * t55) * t78;
t13 = (-t43 * t51 + t44 * t54) * t77;
t12 = Ifges(2,3) * t24 + t15 * t65;
t11 = Ifges(2,3) * t23 + t14 * t64;
t10 = Ifges(2,3) * t22 + t13 * t63;
t9 = t15 * t96 + t24 * t65;
t8 = t14 * t96 + t23 * t64;
t7 = t13 * t96 + t22 * t63;
t6 = (mrSges(2,1) * t62 + mrSges(2,2) * t59) * t89 - t92 * (mrSges(2,1) * t59 - mrSges(2,2) * t62) + (-Ifges(2,3) * t92 + pkin(2) * t65) * t79 * t18;
t5 = (mrSges(2,1) * t61 + mrSges(2,2) * t58) * t88 - t91 * (mrSges(2,1) * t58 - mrSges(2,2) * t61) + (-Ifges(2,3) * t91 + pkin(2) * t64) * t78 * t17;
t4 = (mrSges(2,1) * t60 + mrSges(2,2) * t57) * t87 - t90 * (mrSges(2,1) * t57 - mrSges(2,2) * t60) + (-Ifges(2,3) * t90 + pkin(2) * t63) * t77 * t16;
t3 = -t58 * t96 + (-t88 * mrSges(2,1) - t91 * mrSges(2,2) + (-t64 * t91 + t116) * t78) * t17;
t2 = -t57 * t96 + (-t87 * mrSges(2,1) - t90 * mrSges(2,2) + (-t63 * t90 + t116) * t77) * t16;
t1 = -t59 * t96 + (-t89 * mrSges(2,1) - t92 * mrSges(2,2) + (-t65 * t92 + t116) * t79) * t18;
t72 = [-t49 * t84 - t80 * t50 + (t86 - g(1)) * m(3) + ((-t35 * t110 + t29 * t47) * t42 + (-t35 * t113 + t29 * t48) * t39 + t47 * t1 - t6 * t110) * t79 + ((-t33 * t111 + t27 * t45) * t41 + (-t33 * t114 + t27 * t46) * t38 + t45 * t3 - t5 * t111) * t78 + ((-t31 * t112 + t25 * t43) * t40 + (-t31 * t115 + t25 * t44) * t37 + t43 * t2 - t4 * t112) * t77; -t80 * t49 + t50 * t84 + (t85 - g(2)) * m(3) + ((-t36 * t110 + t30 * t47) * t42 + (-t36 * t113 + t30 * t48) * t39 + t48 * t1 - t6 * t113) * t79 + ((-t34 * t111 + t28 * t45) * t41 + (-t34 * t114 + t28 * t46) * t38 + t46 * t3 - t5 * t114) * t78 + ((-t32 * t112 + t26 * t43) * t40 + (-t32 * t115 + t26 * t44) * t37 + t44 * t2 - t4 * t115) * t77; t15 * t1 + t24 * t6 + t14 * t3 + t23 * t5 + t13 * t2 + t22 * t4 - t49 * t86 + t50 * t85 + Ifges(3,3) * t84 - (-g(1) * t99 - g(2) * t98) * t75 + t76 * (g(1) * t98 - g(2) * t99) + ((-t12 * t110 + t47 * t9) * t42 + (-t12 * t113 + t48 * t9) * t39) * t79 + ((-t11 * t111 + t45 * t8) * t41 + (-t11 * t114 + t46 * t8) * t38) * t78 + ((-t10 * t112 + t43 * t7) * t40 + (-t10 * t115 + t44 * t7) * t37) * t77;];
tauX  = t72;
