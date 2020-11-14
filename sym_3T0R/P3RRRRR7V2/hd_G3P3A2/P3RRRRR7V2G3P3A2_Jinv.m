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
% Datum: 2020-08-07 10:55
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G3P3A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3P3A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3P3A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3P3A2_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3P3A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3P3A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:53:08
% EndTime: 2020-08-07 10:53:08
% DurationCPUTime: 0.42s
% Computational Cost: add. (267->157), mult. (429->183), div. (21->10), fcn. (306->63), ass. (0->105)
t132 = 2 * pkin(3);
t43 = (pkin(5) + pkin(6) + pkin(7));
t131 = 2 * t43;
t66 = cos(qJ(2,3));
t130 = 0.2e1 * t66 ^ 2;
t69 = cos(qJ(2,2));
t129 = 0.2e1 * t69 ^ 2;
t72 = cos(qJ(2,1));
t128 = 0.2e1 * t72 ^ 2;
t65 = cos(qJ(3,3));
t47 = t65 ^ 2;
t127 = pkin(3) * t47;
t68 = cos(qJ(3,2));
t49 = t68 ^ 2;
t126 = pkin(3) * t49;
t71 = cos(qJ(3,1));
t51 = t71 ^ 2;
t125 = pkin(3) * t51;
t56 = sin(qJ(3,3));
t124 = t56 * pkin(1);
t59 = sin(qJ(3,2));
t123 = t59 * pkin(1);
t62 = sin(qJ(3,1));
t122 = t62 * pkin(1);
t121 = t65 * pkin(2);
t120 = t68 * pkin(2);
t119 = t71 * pkin(2);
t67 = cos(qJ(1,3));
t75 = -pkin(3) / 0.2e1;
t118 = (t127 + t121 / 0.2e1 + t75) * t67;
t70 = cos(qJ(1,2));
t117 = (t126 + t120 / 0.2e1 + t75) * t70;
t73 = cos(qJ(1,1));
t116 = (t125 + t119 / 0.2e1 + t75) * t73;
t37 = t65 * pkin(3);
t76 = pkin(2) / 0.2e1;
t115 = (t37 + t76) * t56;
t38 = t68 * pkin(3);
t114 = (t38 + t76) * t59;
t39 = t71 * pkin(3);
t113 = (t39 + t76) * t62;
t57 = sin(qJ(2,3));
t112 = t56 * t57;
t60 = sin(qJ(2,2));
t111 = t59 * t60;
t63 = sin(qJ(2,1));
t110 = t62 * t63;
t109 = t65 * (pkin(1) * t57 - pkin(3) * t56);
t108 = t68 * (pkin(1) * t60 - pkin(3) * t59);
t107 = t71 * (pkin(1) * t63 - pkin(3) * t62);
t58 = sin(qJ(1,3));
t106 = pkin(1) * t67 + t58 * t43;
t61 = sin(qJ(1,2));
t105 = pkin(1) * t70 + t61 * t43;
t64 = sin(qJ(1,1));
t104 = pkin(1) * t73 + t64 * t43;
t103 = -qJ(3,1) + qJ(1,1);
t102 = qJ(3,1) + qJ(1,1);
t101 = -qJ(3,2) + qJ(1,2);
t100 = qJ(3,2) + qJ(1,2);
t99 = -qJ(3,3) + qJ(1,3);
t98 = qJ(3,3) + qJ(1,3);
t97 = 0.2e1 * pkin(1);
t87 = 0.1e1 / pkin(2);
t96 = 0.1e1 / ((t37 + pkin(2)) * t66 - pkin(3) * t112 + pkin(1)) / t56 * t87;
t95 = 0.1e1 / ((t38 + pkin(2)) * t69 - pkin(3) * t111 + pkin(1)) / t59 * t87;
t94 = 0.1e1 / ((t39 + pkin(2)) * t72 - pkin(3) * t110 + pkin(1)) / t62 * t87;
t93 = t67 * t112;
t92 = t70 * t111;
t91 = t73 * t110;
t90 = pkin(2) * t93 + (t132 * t93 - t106) * t65;
t89 = pkin(2) * t92 + (t132 * t92 - t105) * t68;
t88 = pkin(2) * t91 + (t132 * t91 - t104) * t71;
t86 = pkin(2) ^ 2;
t85 = -0.2e1 * qJ(2,1);
t84 = 0.2e1 * qJ(2,1);
t83 = 0.2e1 * qJ(3,1);
t82 = -0.2e1 * qJ(2,2);
t81 = 0.2e1 * qJ(2,2);
t80 = 0.2e1 * qJ(3,2);
t79 = -0.2e1 * qJ(2,3);
t78 = 0.2e1 * qJ(2,3);
t77 = 0.2e1 * qJ(3,3);
t55 = legFrame(1,2);
t54 = legFrame(2,2);
t53 = legFrame(3,2);
t36 = cos(t55);
t35 = cos(t54);
t34 = cos(t53);
t33 = sin(t55);
t32 = sin(t54);
t31 = sin(t53);
t30 = -qJ(2,1) + t103;
t29 = qJ(2,1) + t102;
t28 = -qJ(2,2) + t101;
t27 = qJ(2,2) + t100;
t26 = -qJ(2,3) + t99;
t25 = qJ(2,3) + t98;
t6 = t122 + (-pkin(3) + t119 + 0.2e1 * t125) * t63;
t5 = t123 + (-pkin(3) + t120 + 0.2e1 * t126) * t60;
t4 = t124 + (-pkin(3) + t121 + 0.2e1 * t127) * t57;
t3 = t104 * t110 + (t51 - 0.1e1) * t73 * pkin(3);
t2 = t105 * t111 + (t49 - 0.1e1) * t70 * pkin(3);
t1 = t106 * t112 + (t47 - 0.1e1) * t67 * pkin(3);
t7 = [((t113 * t33 + t116 * t36) * t128 + (t33 * t6 - t36 * t88) * t72 - t3 * t36 + t33 * t107) * t94, ((t113 * t36 - t116 * t33) * t128 + (t33 * t88 + t36 * t6) * t72 + t3 * t33 + t36 * t107) * t94, ((-sin(t30) - sin(t29)) * t97 + (cos(t30) + cos(t29)) * t131 + (-sin(qJ(1,1) + t85 - 0.2e1 * qJ(3,1)) - sin(qJ(1,1) + t84 + t83) - 0.2e1 * t64) * pkin(3) + (-sin(t85 + t103) - sin(t84 + t102) - sin(t103) - sin(t102)) * pkin(2)) / (-t86 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (pkin(2) * sin(qJ(2,1) + qJ(3,1)) + 0.2e1 * t122 + (sin(t83 + qJ(2,1)) - t63) * pkin(3))) / 0.2e1; ((t114 * t32 + t117 * t35) * t129 + (t32 * t5 - t35 * t89) * t69 - t2 * t35 + t32 * t108) * t95, ((t114 * t35 - t117 * t32) * t129 + (t32 * t89 + t35 * t5) * t69 + t2 * t32 + t35 * t108) * t95, ((-sin(t28) - sin(t27)) * t97 + (cos(t28) + cos(t27)) * t131 + (-sin(qJ(1,2) + t82 - 0.2e1 * qJ(3,2)) - sin(qJ(1,2) + t81 + t80) - 0.2e1 * t61) * pkin(3) + (-sin(t82 + t101) - sin(t81 + t100) - sin(t101) - sin(t100)) * pkin(2)) / (-t86 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (pkin(2) * sin(qJ(2,2) + qJ(3,2)) + 0.2e1 * t123 + (sin(t80 + qJ(2,2)) - t60) * pkin(3))) / 0.2e1; ((t115 * t31 + t118 * t34) * t130 + (t31 * t4 - t34 * t90) * t66 - t1 * t34 + t31 * t109) * t96, ((t115 * t34 - t118 * t31) * t130 + (t31 * t90 + t34 * t4) * t66 + t1 * t31 + t34 * t109) * t96, ((-sin(t25) - sin(t26)) * t97 + (cos(t25) + cos(t26)) * t131 + (-sin(qJ(1,3) + t78 + t77) - sin(qJ(1,3) + t79 - 0.2e1 * qJ(3,3)) - 0.2e1 * t58) * pkin(3) + (-sin(t79 + t99) - sin(t78 + t98) - sin(t98) - sin(t99)) * pkin(2)) / (-t86 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (pkin(2) * sin(qJ(2,3) + qJ(3,3)) + 0.2e1 * t124 + (sin(t77 + qJ(2,3)) - t57) * pkin(3))) / 0.2e1;];
Jinv  = t7;
