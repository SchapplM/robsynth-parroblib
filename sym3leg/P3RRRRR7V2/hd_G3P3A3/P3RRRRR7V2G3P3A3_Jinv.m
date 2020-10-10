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
% Datum: 2020-08-07 10:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G3P3A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3P3A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3P3A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3P3A3_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3P3A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3P3A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:55:40
% EndTime: 2020-08-07 10:55:41
% DurationCPUTime: 0.34s
% Computational Cost: add. (222->90), mult. (480->168), div. (36->8), fcn. (306->24), ass. (0->103)
t60 = cos(qJ(2,3));
t124 = 0.2e1 * t60 ^ 2;
t63 = cos(qJ(2,2));
t123 = 0.2e1 * t63 ^ 2;
t66 = cos(qJ(2,1));
t122 = 0.2e1 * t66 ^ 2;
t59 = cos(qJ(3,3));
t39 = t59 ^ 2;
t61 = cos(qJ(1,3));
t35 = pkin(5) + pkin(6) + pkin(7);
t52 = sin(qJ(1,3));
t90 = pkin(1) * t61 + t52 * t35;
t50 = sin(qJ(3,3));
t51 = sin(qJ(2,3));
t97 = t50 * t51;
t121 = (t90 * t97 + (t39 - 0.1e1) * t61 * pkin(3)) * pkin(3);
t62 = cos(qJ(3,2));
t41 = t62 ^ 2;
t64 = cos(qJ(1,2));
t55 = sin(qJ(1,2));
t89 = pkin(1) * t64 + t55 * t35;
t53 = sin(qJ(3,2));
t54 = sin(qJ(2,2));
t96 = t53 * t54;
t120 = (t89 * t96 + (t41 - 0.1e1) * t64 * pkin(3)) * pkin(3);
t65 = cos(qJ(3,1));
t43 = t65 ^ 2;
t67 = cos(qJ(1,1));
t58 = sin(qJ(1,1));
t88 = pkin(1) * t67 + t58 * t35;
t56 = sin(qJ(3,1));
t57 = sin(qJ(2,1));
t95 = t56 * t57;
t119 = (t88 * t95 + (t43 - 0.1e1) * t67 * pkin(3)) * pkin(3);
t118 = t50 * pkin(3);
t117 = t53 * pkin(3);
t116 = t56 * pkin(3);
t115 = t59 * pkin(3);
t114 = t62 * pkin(3);
t113 = t65 * pkin(3);
t68 = pkin(3) ^ 2;
t100 = t39 * t68;
t86 = pkin(2) * t115;
t70 = pkin(2) ^ 2;
t87 = -t68 / 0.2e1 + t70 / 0.2e1;
t13 = t86 + t87 + t100;
t112 = t13 * t61;
t85 = pkin(2) * t114;
t99 = t41 * t68;
t14 = t85 + t87 + t99;
t111 = t14 * t64;
t84 = pkin(2) * t113;
t98 = t43 * t68;
t15 = t84 + t87 + t98;
t110 = t15 * t67;
t22 = pkin(2) + t115;
t47 = legFrame(3,2);
t25 = sin(t47);
t109 = t22 * t25;
t28 = cos(t47);
t108 = t22 * t28;
t107 = t22 * t60;
t23 = pkin(2) + t114;
t48 = legFrame(2,2);
t26 = sin(t48);
t106 = t23 * t26;
t29 = cos(t48);
t105 = t23 * t29;
t104 = t23 * t63;
t24 = pkin(2) + t113;
t49 = legFrame(1,2);
t27 = sin(t49);
t103 = t24 * t27;
t30 = cos(t49);
t102 = t24 * t30;
t101 = t24 * t66;
t94 = t61 * t35;
t93 = t64 * t35;
t92 = t67 * t35;
t91 = 0.1e1 / pkin(3) / pkin(2);
t83 = t22 * t118;
t82 = t23 * t117;
t81 = t24 * t116;
t80 = pkin(3) * t97;
t79 = pkin(3) * t96;
t78 = pkin(3) * t95;
t77 = -0.2e1 * t80;
t76 = -0.2e1 * t79;
t75 = -0.2e1 * t78;
t74 = 0.1e1 / (pkin(1) - t80 + t107) / t50 * t91;
t73 = 0.1e1 / (pkin(1) - t79 + t104) / t53 * t91;
t72 = 0.1e1 / (pkin(1) - t78 + t101) / t56 * t91;
t34 = -t68 + t70;
t18 = pkin(1) * t57 - t116;
t17 = pkin(1) * t54 - t117;
t16 = pkin(1) * t51 - t118;
t9 = t67 * t75 + t88;
t8 = t64 * t76 + t89;
t7 = t61 * t77 + t90;
t6 = pkin(1) * t116 + (t34 + 0.2e1 * t84 + 0.2e1 * t98) * t57;
t5 = pkin(1) * t117 + (t34 + 0.2e1 * t85 + 0.2e1 * t99) * t54;
t4 = pkin(1) * t118 + (t34 + 0.2e1 * t86 + 0.2e1 * t100) * t51;
t1 = [((-t110 * t30 - t27 * t81) * t122 + (-t102 * t9 - t6 * t27) * t66 + t30 * t119 - t18 * t103) * t72, ((t110 * t27 - t30 * t81) * t122 + (t103 * t9 - t6 * t30) * t66 - t27 * t119 - t18 * t102) * t72, (t15 * t58 * t122 + ((pkin(1) + t75) * t58 - t92) * t101 - ((pkin(1) * t95 + pkin(3) * t43 - pkin(3)) * t58 - t92 * t95) * pkin(3)) * t72; ((-t111 * t29 - t26 * t82) * t123 + (-t105 * t8 - t5 * t26) * t63 + t29 * t120 - t17 * t106) * t73, ((t111 * t26 - t29 * t82) * t123 + (t106 * t8 - t5 * t29) * t63 - t26 * t120 - t17 * t105) * t73, (t14 * t55 * t123 + ((pkin(1) + t76) * t55 - t93) * t104 - ((pkin(1) * t96 + pkin(3) * t41 - pkin(3)) * t55 - t93 * t96) * pkin(3)) * t73; ((-t112 * t28 - t25 * t83) * t124 + (-t108 * t7 - t4 * t25) * t60 + t28 * t121 - t16 * t109) * t74, ((t112 * t25 - t28 * t83) * t124 + (t109 * t7 - t4 * t28) * t60 - t25 * t121 - t16 * t108) * t74, (t13 * t52 * t124 + ((pkin(1) + t77) * t52 - t94) * t107 - ((pkin(1) * t97 + pkin(3) * t39 - pkin(3)) * t52 - t94 * t97) * pkin(3)) * t74;];
Jinv  = t1;
