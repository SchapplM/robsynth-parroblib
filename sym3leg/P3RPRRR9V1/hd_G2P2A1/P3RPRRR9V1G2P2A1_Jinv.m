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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:56
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR9V1G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2P2A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:55:58
% EndTime: 2020-08-06 18:55:58
% DurationCPUTime: 0.20s
% Computational Cost: add. (174->83), mult. (321->141), div. (15->6), fcn. (267->23), ass. (0->76)
t94 = 2 * pkin(3);
t44 = cos(pkin(7));
t93 = 0.2e1 * t44 ^ 2;
t92 = pkin(5) + pkin(6);
t54 = cos(qJ(3,3));
t40 = t54 ^ 2;
t91 = t40 * pkin(3);
t56 = cos(qJ(3,2));
t41 = t56 ^ 2;
t90 = t41 * pkin(3);
t58 = cos(qJ(3,1));
t42 = t58 ^ 2;
t89 = t42 * pkin(3);
t88 = t54 * pkin(2);
t87 = t56 * pkin(2);
t86 = t58 * pkin(2);
t36 = qJ(2,3) + t92;
t24 = 0.1e1 / t36;
t43 = sin(pkin(7));
t48 = sin(qJ(3,3));
t76 = t43 * t48;
t85 = 0.1e1 / (t44 * t54 - t76) * t24;
t37 = qJ(2,2) + t92;
t25 = 0.1e1 / t37;
t50 = sin(qJ(3,2));
t75 = t43 * t50;
t84 = 0.1e1 / (t44 * t56 - t75) * t25;
t38 = qJ(2,1) + t92;
t26 = 0.1e1 / t38;
t52 = sin(qJ(3,1));
t74 = t43 * t52;
t83 = 0.1e1 / (t44 * t58 - t74) * t26;
t49 = sin(qJ(1,3));
t60 = -pkin(3) / 0.2e1;
t82 = (t91 + t88 / 0.2e1 + t60) * t49;
t51 = sin(qJ(1,2));
t81 = (t90 + t87 / 0.2e1 + t60) * t51;
t53 = sin(qJ(1,1));
t80 = (t89 + t86 / 0.2e1 + t60) * t53;
t61 = pkin(2) / 0.2e1;
t79 = (t54 * pkin(3) + t61) * t48;
t78 = (t56 * pkin(3) + t61) * t50;
t77 = (t58 * pkin(3) + t61) * t52;
t23 = pkin(1) * t43;
t73 = t54 * (-t48 * pkin(3) + t23);
t72 = t56 * (-t50 * pkin(3) + t23);
t71 = t58 * (-t52 * pkin(3) + t23);
t70 = t49 * t76;
t69 = t51 * t75;
t68 = t53 * t74;
t55 = cos(qJ(1,3));
t67 = pkin(1) * t49 - t55 * t36;
t57 = cos(qJ(1,2));
t66 = pkin(1) * t51 - t57 * t37;
t59 = cos(qJ(1,1));
t65 = pkin(1) * t53 - t59 * t38;
t64 = pkin(2) * t70 + (t70 * t94 - t67) * t54;
t63 = pkin(2) * t69 + (t69 * t94 - t66) * t56;
t62 = pkin(2) * t68 + (t68 * t94 - t65) * t58;
t47 = legFrame(1,2);
t46 = legFrame(2,2);
t45 = legFrame(3,2);
t32 = cos(t47);
t31 = cos(t46);
t30 = cos(t45);
t29 = sin(t47);
t28 = sin(t46);
t27 = sin(t45);
t22 = t44 * pkin(2) + pkin(1);
t6 = pkin(1) * t52 + (-pkin(3) + t86 + 0.2e1 * t89) * t43;
t5 = pkin(1) * t50 + (-pkin(3) + t87 + 0.2e1 * t90) * t43;
t4 = pkin(1) * t48 + (-pkin(3) + t88 + 0.2e1 * t91) * t43;
t3 = t65 * t74 + (t42 - 0.1e1) * t53 * pkin(3);
t2 = t66 * t75 + (t41 - 0.1e1) * t51 * pkin(3);
t1 = t67 * t76 + (t40 - 0.1e1) * t49 * pkin(3);
t7 = [((t29 * t77 + t32 * t80) * t93 + (t29 * t6 - t62 * t32) * t44 - t3 * t32 + t29 * t71) * t83, ((-t29 * t80 + t32 * t77) * t93 + (t62 * t29 + t32 * t6) * t44 + t3 * t29 + t32 * t71) * t83, (t53 * t38 + (pkin(3) * cos(pkin(7) + qJ(3,1)) + t22) * t59) * t26; ((t28 * t78 + t31 * t81) * t93 + (t28 * t5 - t63 * t31) * t44 - t2 * t31 + t28 * t72) * t84, ((-t28 * t81 + t31 * t78) * t93 + (t63 * t28 + t31 * t5) * t44 + t2 * t28 + t31 * t72) * t84, (t51 * t37 + (cos(pkin(7) + qJ(3,2)) * pkin(3) + t22) * t57) * t25; ((t27 * t79 + t30 * t82) * t93 + (t27 * t4 - t64 * t30) * t44 - t1 * t30 + t27 * t73) * t85, ((-t27 * t82 + t30 * t79) * t93 + (t64 * t27 + t30 * t4) * t44 + t1 * t27 + t30 * t73) * t85, (t49 * t36 + (cos(pkin(7) + qJ(3,3)) * pkin(3) + t22) * t55) * t24;];
Jinv  = t7;
