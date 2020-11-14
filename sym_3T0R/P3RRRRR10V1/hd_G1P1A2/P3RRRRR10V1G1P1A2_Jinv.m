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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
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
% Datum: 2020-08-06 22:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V1G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G1P1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G1P1A2_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 22:13:11
% EndTime: 2020-08-06 22:13:12
% DurationCPUTime: 0.62s
% Computational Cost: add. (246->127), mult. (648->228), div. (18->9), fcn. (633->26), ass. (0->102)
t44 = sin(pkin(3));
t45 = cos(pkin(3));
t93 = t44 * t45;
t113 = pkin(1) * t45;
t58 = cos(qJ(3,3));
t112 = pkin(2) * t58 ^ 2;
t61 = cos(qJ(3,2));
t111 = pkin(2) * t61 ^ 2;
t64 = cos(qJ(3,1));
t110 = pkin(2) * t64 ^ 2;
t49 = sin(qJ(3,3));
t109 = pkin(2) * t49;
t52 = sin(qJ(3,2));
t108 = pkin(2) * t52;
t55 = sin(qJ(3,1));
t107 = pkin(2) * t55;
t106 = pkin(2) * t58;
t105 = pkin(2) * t61;
t104 = pkin(2) * t64;
t103 = pkin(5) * t44;
t19 = pkin(1) * t109;
t36 = 0.1e1 / t58;
t50 = sin(qJ(2,3));
t59 = cos(qJ(2,3));
t89 = t58 * t59;
t97 = t50 * pkin(6);
t102 = 0.1e1 / ((t59 * pkin(6) - t50 * t106) * t113 + (t19 + (-pkin(2) * t89 - t97) * pkin(5)) * t44) * t36;
t20 = pkin(1) * t108;
t39 = 0.1e1 / t61;
t53 = sin(qJ(2,2));
t62 = cos(qJ(2,2));
t88 = t61 * t62;
t96 = t53 * pkin(6);
t101 = 0.1e1 / ((t62 * pkin(6) - t53 * t105) * t113 + (t20 + (-pkin(2) * t88 - t96) * pkin(5)) * t44) * t39;
t21 = pkin(1) * t107;
t42 = 0.1e1 / t64;
t56 = sin(qJ(2,1));
t65 = cos(qJ(2,1));
t87 = t64 * t65;
t95 = t56 * pkin(6);
t100 = 0.1e1 / ((t65 * pkin(6) - t56 * t104) * t113 + (t21 + (-pkin(2) * t87 - t95) * pkin(5)) * t44) * t42;
t34 = t45 ^ 2;
t16 = t49 * pkin(5) + pkin(2);
t4 = pkin(6) * t89 + (t16 - 0.2e1 * t112) * t50;
t99 = t4 * t34;
t17 = t52 * pkin(5) + pkin(2);
t5 = pkin(6) * t88 + (t17 - 0.2e1 * t111) * t53;
t98 = t5 * t34;
t18 = t55 * pkin(5) + pkin(2);
t6 = pkin(6) * t87 + (t18 - 0.2e1 * t110) * t56;
t94 = t6 * t34;
t92 = t44 * t49;
t91 = t44 * t52;
t90 = t44 * t55;
t86 = pkin(6) * t113;
t85 = t50 * t112;
t84 = t53 * t111;
t83 = t56 * t110;
t46 = legFrame(3,3);
t22 = sin(t46);
t25 = cos(t46);
t51 = sin(qJ(1,3));
t60 = cos(qJ(1,3));
t10 = t22 * t60 + t51 * t25;
t82 = t10 * t92;
t47 = legFrame(2,3);
t23 = sin(t47);
t26 = cos(t47);
t54 = sin(qJ(1,2));
t63 = cos(qJ(1,2));
t11 = t23 * t63 + t54 * t26;
t81 = t11 * t91;
t48 = legFrame(1,3);
t24 = sin(t48);
t27 = cos(t48);
t57 = sin(qJ(1,1));
t66 = cos(qJ(1,1));
t12 = t24 * t66 + t57 * t27;
t80 = t12 * t90;
t79 = t49 * t59 * t50;
t78 = t52 * t62 * t53;
t77 = t55 * t65 * t56;
t76 = pkin(6) * (t34 - 0.1e1);
t75 = pkin(1) + t97;
t74 = pkin(1) + t96;
t73 = pkin(1) + t95;
t72 = pkin(6) * (t59 - 0.1e1) * (t59 + 0.1e1) * t92;
t71 = pkin(6) * (t62 - 0.1e1) * (t62 + 0.1e1) * t91;
t70 = pkin(6) * (t65 - 0.1e1) * (t65 + 0.1e1) * t90;
t69 = t44 * t79 * t106;
t68 = t44 * t78 * t105;
t67 = t44 * t77 * t104;
t43 = t65 ^ 2;
t40 = t62 ^ 2;
t37 = t59 ^ 2;
t15 = (t43 - 0.2e1) * t107 - pkin(5);
t14 = (t40 - 0.2e1) * t108 - pkin(5);
t13 = (t37 - 0.2e1) * t109 - pkin(5);
t9 = t24 * t57 - t66 * t27;
t8 = t23 * t54 - t63 * t26;
t7 = t22 * t51 - t60 * t25;
t1 = [(t12 * t94 + ((-t12 * t15 * t44 - t73 * t9) * t64 + (-t9 * t110 - t80 * t95) * t65) * t45 + t12 * t83 + t9 * t67 + (pkin(1) * t9 * t90 - t18 * t12) * t56 - t9 * t70) * t100, (-pkin(1) * t80 * t56 + (t94 + (-pkin(6) * t77 - t15 * t64) * t93 + t83 - t56 * t18) * t9 + ((t65 * t110 + t73 * t64) * t45 - t67 + t70) * t12) * t100, ((t15 * t34 + pkin(5)) * t64 + ((-t43 + 0.1e1) * t104 + (t56 * t76 - pkin(1)) * t65) * t55 + t6 * t93) * t42 / ((t65 * t103 + t56 * t113) * t104 - t65 * t86 - t44 * (-pkin(5) * t95 + t21)); (t11 * t98 + ((-t11 * t14 * t44 - t74 * t8) * t61 + (-t8 * t111 - t81 * t96) * t62) * t45 + t11 * t84 + t8 * t68 + (pkin(1) * t8 * t91 - t17 * t11) * t53 - t8 * t71) * t101, (-pkin(1) * t81 * t53 + (t98 + (-pkin(6) * t78 - t14 * t61) * t93 + t84 - t53 * t17) * t8 + ((t62 * t111 + t74 * t61) * t45 - t68 + t71) * t11) * t101, ((t14 * t34 + pkin(5)) * t61 + ((-t40 + 0.1e1) * t105 + (t53 * t76 - pkin(1)) * t62) * t52 + t5 * t93) * t39 / ((t62 * t103 + t53 * t113) * t105 - t62 * t86 - t44 * (-pkin(5) * t96 + t20)); (t10 * t99 + ((-t10 * t13 * t44 - t75 * t7) * t58 + (-t7 * t112 - t82 * t97) * t59) * t45 + t10 * t85 + t7 * t69 + (pkin(1) * t7 * t92 - t16 * t10) * t50 - t7 * t72) * t102, (-pkin(1) * t82 * t50 + (t99 + (-pkin(6) * t79 - t13 * t58) * t93 + t85 - t50 * t16) * t7 + ((t59 * t112 + t75 * t58) * t45 - t69 + t72) * t10) * t102, ((t13 * t34 + pkin(5)) * t58 + ((-t37 + 0.1e1) * t106 + (t50 * t76 - pkin(1)) * t59) * t49 + t4 * t93) * t36 / ((t59 * t103 + t50 * t113) * t106 - t59 * t86 - t44 * (-pkin(5) * t97 + t19));];
Jinv  = t1;
