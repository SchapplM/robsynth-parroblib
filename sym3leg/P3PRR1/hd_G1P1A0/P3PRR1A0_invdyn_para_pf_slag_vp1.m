% Calculate vector of inverse dynamics forces for parallel robot
% P3PRR1G1P1A0
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
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
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

function tauX = P3PRR1G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:34
% EndTime: 2019-05-03 14:47:34
% DurationCPUTime: 0.65s
% Computational Cost: add. (845->166), mult. (1565->320), div. (450->4), fcn. (1400->14), ass. (0->125)
t85 = sin(qJ(2,3));
t88 = cos(qJ(2,3));
t63 = rSges(2,1) * t88 - rSges(2,2) * t85;
t124 = m(2) * t63;
t86 = sin(qJ(2,2));
t89 = cos(qJ(2,2));
t64 = rSges(2,1) * t89 - rSges(2,2) * t86;
t123 = m(2) * t64;
t87 = sin(qJ(2,1));
t90 = cos(qJ(2,1));
t65 = rSges(2,1) * t90 - rSges(2,2) * t87;
t122 = m(2) * t65;
t94 = m(1) + m(2);
t121 = pkin(2) * t94;
t104 = 0.1e1 / pkin(2);
t120 = m(2) * t104;
t75 = 0.1e1 / t85;
t119 = t63 * t75;
t76 = 0.1e1 / t86;
t118 = t64 * t76;
t77 = 0.1e1 / t87;
t117 = t65 * t77;
t79 = legFrame(3,3);
t67 = sin(t79);
t116 = t104 * t67;
t80 = legFrame(2,3);
t68 = sin(t80);
t115 = t104 * t68;
t81 = legFrame(1,3);
t69 = sin(t81);
t114 = t104 * t69;
t70 = cos(t79);
t113 = t104 * t70;
t71 = cos(t80);
t112 = t104 * t71;
t72 = cos(t81);
t111 = t104 * t72;
t110 = t104 * t75;
t109 = t104 * t76;
t108 = t104 * t77;
t107 = t63 * t120;
t106 = t64 * t120;
t105 = t65 * t120;
t103 = koppelP(1,1);
t102 = koppelP(2,1);
t101 = koppelP(3,1);
t100 = koppelP(1,2);
t99 = koppelP(2,2);
t98 = koppelP(3,2);
t97 = rSges(3,1);
t96 = rSges(3,2);
t95 = xP(3);
t93 = xDP(1);
t92 = xDP(2);
t91 = xDP(3);
t84 = xDDP(1);
t83 = xDDP(2);
t82 = xDDP(3);
t78 = t91 ^ 2;
t74 = cos(t95);
t73 = sin(t95);
t66 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3);
t62 = g(1) * t72 + g(2) * t69;
t61 = g(1) * t71 + g(2) * t68;
t60 = g(1) * t70 + g(2) * t67;
t59 = -g(1) * t69 + g(2) * t72;
t58 = -g(1) * t68 + g(2) * t71;
t57 = -g(1) * t67 + g(2) * t70;
t56 = -t100 * t73 + t103 * t74;
t55 = t102 * t74 - t73 * t99;
t54 = t101 * t74 - t73 * t98;
t53 = t100 * t74 + t103 * t73;
t52 = t102 * t73 + t74 * t99;
t51 = t101 * t73 + t74 * t98;
t50 = -t73 * t96 + t74 * t97;
t49 = t73 * t97 + t74 * t96;
t48 = t69 * t90 + t72 * t87;
t47 = -t69 * t87 + t72 * t90;
t46 = t68 * t89 + t71 * t86;
t45 = -t68 * t86 + t71 * t89;
t44 = t67 * t88 + t70 * t85;
t43 = -t67 * t85 + t70 * t88;
t42 = -t53 * t82 - t56 * t78 + t84;
t41 = -t52 * t82 - t55 * t78 + t84;
t40 = -t51 * t82 - t54 * t78 + t84;
t39 = -t53 * t78 + t56 * t82 + t83;
t38 = -t52 * t78 + t55 * t82 + t83;
t37 = -t51 * t78 + t54 * t82 + t83;
t36 = (-t105 * t69 + t48 * t94) * t77;
t35 = (-t105 * t72 + t47 * t94) * t77;
t34 = (-t106 * t68 + t46 * t94) * t76;
t33 = (-t106 * t71 + t45 * t94) * t76;
t32 = (-t107 * t67 + t44 * t94) * t75;
t31 = (-t107 * t70 + t43 * t94) * t75;
t30 = (t53 * t72 - t56 * t69) * t108;
t29 = (t52 * t71 - t55 * t68) * t109;
t28 = (t51 * t70 - t54 * t67) * t110;
t27 = (-t114 * t66 + t122 * t48) * t77;
t26 = (-t111 * t66 + t122 * t47) * t77;
t25 = (-t115 * t66 + t123 * t46) * t76;
t24 = (-t112 * t66 + t123 * t45) * t76;
t23 = (-t116 * t66 + t124 * t44) * t75;
t22 = (-t113 * t66 + t124 * t43) * t75;
t21 = (-t72 * (-t53 * t91 + t93) - t69 * (t56 * t91 + t92)) * t108;
t20 = (-t71 * (-t52 * t91 + t93) - t68 * (t55 * t91 + t92)) * t109;
t19 = (-t70 * (-t51 * t91 + t93) - t67 * (t54 * t91 + t92)) * t110;
t18 = t21 ^ 2;
t17 = t20 ^ 2;
t16 = t19 ^ 2;
t15 = (-t47 * t53 + t48 * t56) * t77;
t14 = (-t45 * t52 + t46 * t55) * t76;
t13 = (-t43 * t51 + t44 * t54) * t75;
t12 = t122 * t30 + t15 * t94;
t11 = t123 * t29 + t14 * t94;
t10 = t124 * t28 + t13 * t94;
t9 = t122 * t15 + t30 * t66;
t8 = t123 * t14 + t29 * t66;
t7 = t124 * t13 + t28 * t66;
t6 = -t66 * t18 * t90 * t77 + (pkin(2) * t18 * t117 + t87 * (rSges(2,1) * t62 + rSges(2,2) * t59) - t90 * (rSges(2,1) * t59 - rSges(2,2) * t62)) * m(2);
t5 = -t66 * t17 * t89 * t76 + (pkin(2) * t17 * t118 + t86 * (rSges(2,1) * t61 + rSges(2,2) * t58) - t89 * (rSges(2,1) * t58 - rSges(2,2) * t61)) * m(2);
t4 = -t66 * t16 * t88 * t75 + (pkin(2) * t16 * t119 + t85 * (rSges(2,1) * t60 + rSges(2,2) * t57) - t88 * (rSges(2,1) * t57 - rSges(2,2) * t60)) * m(2);
t3 = -t59 * t94 + (t77 * t121 + (-t87 * rSges(2,1) + (-rSges(2,2) - t117) * t90) * m(2)) * t18;
t2 = -t58 * t94 + (t76 * t121 + (-t86 * rSges(2,1) + (-rSges(2,2) - t118) * t89) * m(2)) * t17;
t1 = -t57 * t94 + (t75 * t121 + (-t85 * rSges(2,1) + (-rSges(2,2) - t119) * t88) * m(2)) * t16;
t125 = [(-t49 * t82 - t78 * t50 - g(1) + t84) * m(3) + ((-t111 * t26 + t35 * t47) * t42 + (-t114 * t26 + t35 * t48) * t39 + t47 * t3 - t6 * t111) * t77 + ((-t112 * t24 + t33 * t45) * t41 + (-t115 * t24 + t33 * t46) * t38 + t45 * t2 - t5 * t112) * t76 + ((-t113 * t22 + t31 * t43) * t40 + (-t116 * t22 + t31 * t44) * t37 + t43 * t1 - t4 * t113) * t75; (-t78 * t49 + t50 * t82 - g(2) + t83) * m(3) + ((-t111 * t27 + t36 * t47) * t42 + (-t114 * t27 + t36 * t48) * t39 + t48 * t3 - t6 * t114) * t77 + ((-t112 * t25 + t34 * t45) * t41 + (-t115 * t25 + t34 * t46) * t38 + t46 * t2 - t5 * t115) * t76 + ((-t113 * t23 + t32 * t43) * t40 + (-t116 * t23 + t32 * t44) * t37 + t44 * t1 - t4 * t116) * t75; Icges(3,3) * t82 + t13 * t1 + t14 * t2 + t15 * t3 + t28 * t4 + t29 * t5 + t30 * t6 + (-t49 * t84 + t50 * t83 + (t96 ^ 2 + t97 ^ 2) * t82 + (g(1) * t97 + g(2) * t96) * t73 + (g(1) * t96 - g(2) * t97) * t74) * m(3) + ((-t111 * t9 + t12 * t47) * t42 + (-t114 * t9 + t12 * t48) * t39) * t77 + ((t11 * t45 - t112 * t8) * t41 + (t11 * t46 - t115 * t8) * t38) * t76 + ((t10 * t43 - t113 * t7) * t40 + (t10 * t44 - t116 * t7) * t37) * t75;];
tauX  = t125;
