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
% Datum: 2020-08-06 19:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR9V1G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3P3A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:00:48
% EndTime: 2020-08-06 19:00:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (174->83), mult. (321->141), div. (15->6), fcn. (267->23), ass. (0->76)
t97 = 2 * pkin(3);
t47 = cos(pkin(7));
t96 = 0.2e1 * t47 ^ 2;
t95 = pkin(5) + pkin(6);
t57 = cos(qJ(3,3));
t43 = t57 ^ 2;
t94 = t43 * pkin(3);
t59 = cos(qJ(3,2));
t44 = t59 ^ 2;
t93 = t44 * pkin(3);
t61 = cos(qJ(3,1));
t45 = t61 ^ 2;
t92 = t45 * pkin(3);
t91 = t57 * pkin(2);
t90 = t59 * pkin(2);
t89 = t61 * pkin(2);
t39 = qJ(2,3) + t95;
t27 = 0.1e1 / t39;
t46 = sin(pkin(7));
t51 = sin(qJ(3,3));
t79 = t46 * t51;
t88 = 0.1e1 / (t47 * t57 - t79) * t27;
t40 = qJ(2,2) + t95;
t28 = 0.1e1 / t40;
t53 = sin(qJ(3,2));
t78 = t46 * t53;
t87 = 0.1e1 / (t47 * t59 - t78) * t28;
t41 = qJ(2,1) + t95;
t29 = 0.1e1 / t41;
t55 = sin(qJ(3,1));
t77 = t46 * t55;
t86 = 0.1e1 / (t47 * t61 - t77) * t29;
t58 = cos(qJ(1,3));
t63 = -pkin(3) / 0.2e1;
t85 = (t94 + t91 / 0.2e1 + t63) * t58;
t60 = cos(qJ(1,2));
t84 = (t93 + t90 / 0.2e1 + t63) * t60;
t62 = cos(qJ(1,1));
t83 = (t92 + t89 / 0.2e1 + t63) * t62;
t64 = pkin(2) / 0.2e1;
t82 = (pkin(3) * t57 + t64) * t51;
t81 = (pkin(3) * t59 + t64) * t53;
t80 = (pkin(3) * t61 + t64) * t55;
t26 = pkin(1) * t46;
t76 = t57 * (-pkin(3) * t51 + t26);
t75 = t59 * (-pkin(3) * t53 + t26);
t74 = t61 * (-pkin(3) * t55 + t26);
t52 = sin(qJ(1,3));
t73 = pkin(1) * t58 + t52 * t39;
t54 = sin(qJ(1,2));
t72 = pkin(1) * t60 + t54 * t40;
t56 = sin(qJ(1,1));
t71 = pkin(1) * t62 + t56 * t41;
t70 = t58 * t79;
t69 = t60 * t78;
t68 = t62 * t77;
t67 = pkin(2) * t70 + (t70 * t97 - t73) * t57;
t66 = pkin(2) * t69 + (t69 * t97 - t72) * t59;
t65 = pkin(2) * t68 + (t68 * t97 - t71) * t61;
t50 = legFrame(1,2);
t49 = legFrame(2,2);
t48 = legFrame(3,2);
t35 = cos(t50);
t34 = cos(t49);
t33 = cos(t48);
t32 = sin(t50);
t31 = sin(t49);
t30 = sin(t48);
t25 = pkin(2) * t47 + pkin(1);
t6 = pkin(1) * t55 + (-pkin(3) + t89 + 0.2e1 * t92) * t46;
t5 = pkin(1) * t53 + (-pkin(3) + t90 + 0.2e1 * t93) * t46;
t4 = pkin(1) * t51 + (-pkin(3) + t91 + 0.2e1 * t94) * t46;
t3 = t71 * t77 + (t45 - 0.1e1) * t62 * pkin(3);
t2 = t72 * t78 + (t44 - 0.1e1) * t60 * pkin(3);
t1 = t73 * t79 + (t43 - 0.1e1) * t58 * pkin(3);
t7 = [((t32 * t80 + t35 * t83) * t96 + (t32 * t6 - t35 * t65) * t47 - t3 * t35 + t32 * t74) * t86, ((-t32 * t83 + t35 * t80) * t96 + (t32 * t65 + t35 * t6) * t47 + t3 * t32 + t35 * t74) * t86, (t41 * t62 + (-pkin(3) * cos(pkin(7) + qJ(3,1)) - t25) * t56) * t29; ((t31 * t81 + t34 * t84) * t96 + (t31 * t5 - t34 * t66) * t47 - t2 * t34 + t31 * t75) * t87, ((-t31 * t84 + t34 * t81) * t96 + (t31 * t66 + t34 * t5) * t47 + t2 * t31 + t34 * t75) * t87, (t40 * t60 + (-cos(pkin(7) + qJ(3,2)) * pkin(3) - t25) * t54) * t28; ((t30 * t82 + t33 * t85) * t96 + (t30 * t4 - t33 * t67) * t47 - t1 * t33 + t30 * t76) * t88, ((-t30 * t85 + t33 * t82) * t96 + (t30 * t67 + t33 * t4) * t47 + t1 * t30 + t33 * t76) * t88, (t39 * t58 + (-cos(pkin(7) + qJ(3,3)) * pkin(3) - t25) * t52) * t27;];
Jinv  = t7;
