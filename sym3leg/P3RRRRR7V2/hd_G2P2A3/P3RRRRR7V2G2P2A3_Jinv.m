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
% Datum: 2020-08-07 10:19
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G2P2A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2P2A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2P2A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2P2A3_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2P2A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2P2A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:16:41
% EndTime: 2020-08-07 10:16:41
% DurationCPUTime: 0.34s
% Computational Cost: add. (222->90), mult. (480->174), div. (36->8), fcn. (306->24), ass. (0->106)
t57 = cos(qJ(2,3));
t37 = t57 ^ 2;
t121 = 0.2e1 * t37;
t60 = cos(qJ(2,2));
t39 = t60 ^ 2;
t120 = 0.2e1 * t39;
t63 = cos(qJ(2,1));
t41 = t63 ^ 2;
t119 = 0.2e1 * t41;
t56 = cos(qJ(3,3));
t36 = t56 ^ 2;
t49 = sin(qJ(1,3));
t32 = pkin(5) + pkin(6) + pkin(7);
t58 = cos(qJ(1,3));
t74 = pkin(1) * t49 - t58 * t32;
t47 = sin(qJ(3,3));
t48 = sin(qJ(2,3));
t94 = t47 * t48;
t118 = pkin(3) * (t74 * t94 + (t36 - 0.1e1) * t49 * pkin(3));
t59 = cos(qJ(3,2));
t38 = t59 ^ 2;
t52 = sin(qJ(1,2));
t61 = cos(qJ(1,2));
t73 = pkin(1) * t52 - t61 * t32;
t50 = sin(qJ(3,2));
t51 = sin(qJ(2,2));
t92 = t50 * t51;
t117 = pkin(3) * (t73 * t92 + (t38 - 0.1e1) * t52 * pkin(3));
t62 = cos(qJ(3,1));
t40 = t62 ^ 2;
t55 = sin(qJ(1,1));
t64 = cos(qJ(1,1));
t72 = pkin(1) * t55 - t64 * t32;
t53 = sin(qJ(3,1));
t54 = sin(qJ(2,1));
t90 = t53 * t54;
t116 = pkin(3) * (t72 * t90 + (t40 - 0.1e1) * t55 * pkin(3));
t115 = t47 * pkin(3);
t114 = t50 * pkin(3);
t113 = t53 * pkin(3);
t112 = t56 * pkin(3);
t111 = t59 * pkin(3);
t110 = t62 * pkin(3);
t86 = pkin(2) * t112;
t65 = pkin(3) ^ 2;
t67 = pkin(2) ^ 2;
t87 = -t65 / 0.2e1 + t67 / 0.2e1;
t97 = t36 * t65;
t13 = t86 + t87 + t97;
t109 = t13 * t49;
t85 = pkin(2) * t111;
t96 = t38 * t65;
t14 = t85 + t87 + t96;
t108 = t14 * t52;
t84 = pkin(2) * t110;
t95 = t40 * t65;
t15 = t84 + t87 + t95;
t107 = t15 * t55;
t19 = pkin(2) + t112;
t44 = legFrame(3,2);
t22 = sin(t44);
t106 = t19 * t22;
t25 = cos(t44);
t105 = t19 * t25;
t104 = t19 * t57;
t20 = pkin(2) + t111;
t45 = legFrame(2,2);
t23 = sin(t45);
t103 = t20 * t23;
t26 = cos(t45);
t102 = t20 * t26;
t101 = t20 * t60;
t21 = pkin(2) + t110;
t46 = legFrame(1,2);
t24 = sin(t46);
t100 = t21 * t24;
t27 = cos(t46);
t99 = t21 * t27;
t98 = t21 * t63;
t93 = t49 * t32;
t91 = t52 * t32;
t89 = t55 * t32;
t88 = 0.1e1 / pkin(3) / pkin(2);
t83 = t19 * t115;
t82 = t20 * t114;
t81 = t21 * t113;
t80 = pkin(3) * t94;
t79 = pkin(3) * t92;
t78 = pkin(3) * t90;
t77 = -0.2e1 * t80;
t76 = -0.2e1 * t79;
t75 = -0.2e1 * t78;
t71 = 0.1e1 / (pkin(1) - t80 + t104) / t47 * t88;
t70 = 0.1e1 / (pkin(1) - t79 + t101) / t50 * t88;
t69 = 0.1e1 / (pkin(1) - t78 + t98) / t53 * t88;
t31 = -t65 + t67;
t18 = pkin(1) * t54 - t113;
t17 = pkin(1) * t51 - t114;
t16 = pkin(1) * t48 - t115;
t9 = t55 * t75 + t72;
t8 = t52 * t76 + t73;
t7 = t49 * t77 + t74;
t6 = pkin(1) * t113 + (t31 + 0.2e1 * t84 + 0.2e1 * t95) * t54;
t5 = pkin(1) * t114 + (t31 + 0.2e1 * t85 + 0.2e1 * t96) * t51;
t4 = pkin(1) * t115 + (t31 + 0.2e1 * t86 + 0.2e1 * t97) * t48;
t1 = [((-t27 * t107 - t24 * t81) * t119 + (-t24 * t6 - t9 * t99) * t63 + t27 * t116 - t18 * t100) * t69, ((t24 * t107 - t27 * t81) * t119 + (t9 * t100 - t27 * t6) * t63 - t24 * t116 - t18 * t99) * t69, (-0.2e1 * t15 * t64 * t41 - ((pkin(1) + t75) * t64 + t89) * t98 + pkin(3) * ((pkin(1) * t90 + pkin(3) * t40 - pkin(3)) * t64 + t89 * t90)) * t69; ((-t26 * t108 - t23 * t82) * t120 + (-t8 * t102 - t23 * t5) * t60 + t26 * t117 - t17 * t103) * t70, ((t23 * t108 - t26 * t82) * t120 + (t8 * t103 - t26 * t5) * t60 - t23 * t117 - t17 * t102) * t70, (-0.2e1 * t14 * t61 * t39 - ((pkin(1) + t76) * t61 + t91) * t101 + pkin(3) * ((pkin(1) * t92 + pkin(3) * t38 - pkin(3)) * t61 + t91 * t92)) * t70; ((-t25 * t109 - t22 * t83) * t121 + (-t7 * t105 - t22 * t4) * t57 + t25 * t118 - t16 * t106) * t71, ((t22 * t109 - t25 * t83) * t121 + (t7 * t106 - t25 * t4) * t57 - t22 * t118 - t16 * t105) * t71, (-0.2e1 * t13 * t58 * t37 - ((pkin(1) + t77) * t58 + t93) * t104 + pkin(3) * ((pkin(1) * t94 + pkin(3) * t36 - pkin(3)) * t58 + t93 * t94)) * t71;];
Jinv  = t1;
