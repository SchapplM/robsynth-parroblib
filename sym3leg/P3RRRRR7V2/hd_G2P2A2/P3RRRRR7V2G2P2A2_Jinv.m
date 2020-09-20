% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [3x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:16
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G2P2A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2P2A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2P2A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2P2A2_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2P2A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2P2A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:14:18
% EndTime: 2020-08-07 10:14:18
% DurationCPUTime: 0.34s
% Computational Cost: add. (267->157), mult. (429->183), div. (21->10), fcn. (306->63), ass. (0->105)
t129 = 2 * pkin(3);
t40 = (pkin(5) + pkin(6) + pkin(7));
t128 = 2 * t40;
t63 = cos(qJ(2,3));
t127 = 0.2e1 * t63 ^ 2;
t66 = cos(qJ(2,2));
t126 = 0.2e1 * t66 ^ 2;
t69 = cos(qJ(2,1));
t125 = 0.2e1 * t69 ^ 2;
t62 = cos(qJ(3,3));
t44 = t62 ^ 2;
t124 = pkin(3) * t44;
t65 = cos(qJ(3,2));
t46 = t65 ^ 2;
t123 = pkin(3) * t46;
t68 = cos(qJ(3,1));
t48 = t68 ^ 2;
t122 = pkin(3) * t48;
t53 = sin(qJ(3,3));
t121 = t53 * pkin(1);
t56 = sin(qJ(3,2));
t120 = t56 * pkin(1);
t59 = sin(qJ(3,1));
t119 = t59 * pkin(1);
t118 = t62 * pkin(2);
t117 = t65 * pkin(2);
t116 = t68 * pkin(2);
t55 = sin(qJ(1,3));
t72 = -pkin(3) / 0.2e1;
t115 = (t124 + t118 / 0.2e1 + t72) * t55;
t58 = sin(qJ(1,2));
t114 = (t123 + t117 / 0.2e1 + t72) * t58;
t61 = sin(qJ(1,1));
t113 = (t122 + t116 / 0.2e1 + t72) * t61;
t34 = t62 * pkin(3);
t73 = pkin(2) / 0.2e1;
t112 = (t34 + t73) * t53;
t35 = t65 * pkin(3);
t111 = (t35 + t73) * t56;
t36 = t68 * pkin(3);
t110 = (t36 + t73) * t59;
t54 = sin(qJ(2,3));
t109 = t53 * t54;
t57 = sin(qJ(2,2));
t108 = t56 * t57;
t60 = sin(qJ(2,1));
t107 = t59 * t60;
t106 = t62 * (pkin(1) * t54 - t53 * pkin(3));
t105 = t65 * (pkin(1) * t57 - t56 * pkin(3));
t104 = t68 * (pkin(1) * t60 - t59 * pkin(3));
t103 = -qJ(3,1) + qJ(1,1);
t102 = qJ(3,1) + qJ(1,1);
t101 = -qJ(3,2) + qJ(1,2);
t100 = qJ(3,2) + qJ(1,2);
t99 = -qJ(3,3) + qJ(1,3);
t98 = qJ(3,3) + qJ(1,3);
t97 = 0.2e1 * pkin(1);
t84 = 0.1e1 / pkin(2);
t96 = 0.1e1 / ((t34 + pkin(2)) * t63 - pkin(3) * t109 + pkin(1)) / t53 * t84;
t95 = 0.1e1 / ((t35 + pkin(2)) * t66 - pkin(3) * t108 + pkin(1)) / t56 * t84;
t94 = 0.1e1 / ((t36 + pkin(2)) * t69 - pkin(3) * t107 + pkin(1)) / t59 * t84;
t93 = t55 * t109;
t92 = t58 * t108;
t91 = t61 * t107;
t64 = cos(qJ(1,3));
t90 = pkin(1) * t55 - t64 * t40;
t67 = cos(qJ(1,2));
t89 = pkin(1) * t58 - t67 * t40;
t70 = cos(qJ(1,1));
t88 = pkin(1) * t61 - t70 * t40;
t87 = pkin(2) * t93 + (t93 * t129 - t90) * t62;
t86 = pkin(2) * t92 + (t92 * t129 - t89) * t65;
t85 = pkin(2) * t91 + (t91 * t129 - t88) * t68;
t83 = pkin(2) ^ 2;
t82 = -0.2e1 * qJ(2,1);
t81 = 0.2e1 * qJ(2,1);
t80 = 0.2e1 * qJ(3,1);
t79 = -0.2e1 * qJ(2,2);
t78 = 0.2e1 * qJ(2,2);
t77 = 0.2e1 * qJ(3,2);
t76 = -0.2e1 * qJ(2,3);
t75 = 0.2e1 * qJ(2,3);
t74 = 0.2e1 * qJ(3,3);
t52 = legFrame(1,2);
t51 = legFrame(2,2);
t50 = legFrame(3,2);
t33 = cos(t52);
t32 = cos(t51);
t31 = cos(t50);
t30 = sin(t52);
t29 = sin(t51);
t28 = sin(t50);
t27 = -qJ(2,1) + t103;
t26 = qJ(2,1) + t102;
t25 = -qJ(2,2) + t101;
t24 = qJ(2,2) + t100;
t23 = -qJ(2,3) + t99;
t22 = qJ(2,3) + t98;
t6 = t119 + (-pkin(3) + t116 + 0.2e1 * t122) * t60;
t5 = t120 + (-pkin(3) + t117 + 0.2e1 * t123) * t57;
t4 = t121 + (-pkin(3) + t118 + 0.2e1 * t124) * t54;
t3 = t88 * t107 + (t48 - 0.1e1) * t61 * pkin(3);
t2 = t89 * t108 + (t46 - 0.1e1) * t58 * pkin(3);
t1 = t90 * t109 + (t44 - 0.1e1) * t55 * pkin(3);
t7 = [((t30 * t110 + t33 * t113) * t125 + (t30 * t6 - t85 * t33) * t69 - t3 * t33 + t30 * t104) * t94, ((t33 * t110 - t30 * t113) * t125 + (t85 * t30 + t33 * t6) * t69 + t3 * t30 + t33 * t104) * t94, ((cos(t27) + cos(t26)) * t97 + (sin(t27) + sin(t26)) * t128 + (cos(qJ(1,1) + t82 - 0.2e1 * qJ(3,1)) + cos(qJ(1,1) + t81 + t80) + 0.2e1 * t70) * pkin(3) + (cos(t82 + t103) + cos(t81 + t102) + cos(t103) + cos(t102)) * pkin(2)) / (-t83 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (pkin(2) * sin(qJ(2,1) + qJ(3,1)) + 0.2e1 * t119 + (sin(t80 + qJ(2,1)) - t60) * pkin(3))) / 0.2e1; ((t29 * t111 + t32 * t114) * t126 + (t29 * t5 - t86 * t32) * t66 - t2 * t32 + t29 * t105) * t95, ((t32 * t111 - t29 * t114) * t126 + (t86 * t29 + t32 * t5) * t66 + t2 * t29 + t32 * t105) * t95, ((cos(t25) + cos(t24)) * t97 + (sin(t25) + sin(t24)) * t128 + (cos(qJ(1,2) + t79 - 0.2e1 * qJ(3,2)) + cos(qJ(1,2) + t78 + t77) + 0.2e1 * t67) * pkin(3) + (cos(t79 + t101) + cos(t78 + t100) + cos(t101) + cos(t100)) * pkin(2)) / (-t83 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (pkin(2) * sin(qJ(2,2) + qJ(3,2)) + 0.2e1 * t120 + (sin(t77 + qJ(2,2)) - t57) * pkin(3))) / 0.2e1; ((t28 * t112 + t31 * t115) * t127 + (t28 * t4 - t87 * t31) * t63 - t1 * t31 + t28 * t106) * t96, ((t31 * t112 - t28 * t115) * t127 + (t87 * t28 + t31 * t4) * t63 + t1 * t28 + t31 * t106) * t96, ((cos(t22) + cos(t23)) * t97 + (sin(t23) + sin(t22)) * t128 + (cos(qJ(1,3) + t76 - 0.2e1 * qJ(3,3)) + cos(qJ(1,3) + t75 + t74) + 0.2e1 * t64) * pkin(3) + (cos(t76 + t99) + cos(t75 + t98) + cos(t99) + cos(t98)) * pkin(2)) / (-t83 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (pkin(2) * sin(qJ(2,3) + qJ(3,3)) + 0.2e1 * t121 + (sin(t74 + qJ(2,3)) - t54) * pkin(3))) / 0.2e1;];
Jinv  = t7;
