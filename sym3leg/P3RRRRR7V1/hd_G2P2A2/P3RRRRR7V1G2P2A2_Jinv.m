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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-07 08:49
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V1G2P2A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2P2A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2P2A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2P2A2_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2P2A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2P2A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:48:44
% EndTime: 2020-08-07 08:48:44
% DurationCPUTime: 0.28s
% Computational Cost: add. (198->130), mult. (342->153), div. (24->10), fcn. (261->57), ass. (0->90)
t101 = 2 * pkin(2);
t52 = (-pkin(5) - pkin(4));
t100 = 2 * t52;
t62 = 0.1e1 / pkin(1);
t99 = t62 / 0.2e1;
t43 = cos(qJ(3,3));
t98 = t43 * pkin(2);
t46 = cos(qJ(3,2));
t97 = t46 * pkin(2);
t49 = cos(qJ(3,1));
t96 = t49 * pkin(2);
t25 = t43 ^ 2;
t10 = pkin(1) * t43 + t25 * t101 - pkin(2);
t36 = sin(qJ(1,3));
t95 = t10 * t36;
t27 = t46 ^ 2;
t11 = pkin(1) * t46 + t27 * t101 - pkin(2);
t39 = sin(qJ(1,2));
t94 = t11 * t39;
t29 = t49 ^ 2;
t12 = pkin(1) * t49 + t29 * t101 - pkin(2);
t42 = sin(qJ(1,1));
t93 = t12 * t42;
t34 = sin(qJ(3,3));
t92 = (pkin(1) + 0.2e1 * t98) * t34;
t37 = sin(qJ(3,2));
t91 = (pkin(1) + 0.2e1 * t97) * t37;
t40 = sin(qJ(3,1));
t90 = (pkin(1) + 0.2e1 * t96) * t40;
t35 = sin(qJ(2,3));
t89 = t34 * t35;
t88 = t35 * t10;
t38 = sin(qJ(2,2));
t87 = t37 * t38;
t86 = t38 * t11;
t41 = sin(qJ(2,1));
t85 = t40 * t41;
t84 = t41 * t12;
t45 = cos(qJ(1,3));
t83 = t45 * t52;
t48 = cos(qJ(1,2));
t82 = t48 * t52;
t51 = cos(qJ(1,1));
t81 = t51 * t52;
t80 = -qJ(3,1) + qJ(1,1);
t79 = qJ(3,1) + qJ(1,1);
t78 = -qJ(3,2) + qJ(1,2);
t77 = qJ(3,2) + qJ(1,2);
t76 = -qJ(3,3) + qJ(1,3);
t75 = qJ(3,3) + qJ(1,3);
t44 = cos(qJ(2,3));
t74 = 0.1e1 / t34 / ((pkin(1) + t98) * t44 - pkin(2) * t89) * t62;
t47 = cos(qJ(2,2));
t73 = 0.1e1 / t37 / ((pkin(1) + t97) * t47 - pkin(2) * t87) * t62;
t50 = cos(qJ(2,1));
t72 = 0.1e1 / t40 / ((pkin(1) + t96) * t50 - pkin(2) * t85) * t62;
t71 = t34 * t98;
t70 = t37 * t97;
t69 = t40 * t96;
t68 = t36 * t89;
t67 = t39 * t87;
t66 = t42 * t85;
t65 = pkin(1) * t68 + (t68 * t101 - t83) * t43;
t64 = pkin(1) * t67 + (t67 * t101 - t82) * t46;
t63 = pkin(1) * t66 + (t66 * t101 - t81) * t49;
t61 = -0.2e1 * qJ(2,1);
t60 = 0.2e1 * qJ(2,1);
t59 = 0.2e1 * qJ(3,1);
t58 = -0.2e1 * qJ(2,2);
t57 = 0.2e1 * qJ(2,2);
t56 = 0.2e1 * qJ(3,2);
t55 = -0.2e1 * qJ(2,3);
t54 = 0.2e1 * qJ(2,3);
t53 = 0.2e1 * qJ(3,3);
t33 = legFrame(1,2);
t32 = legFrame(2,2);
t31 = legFrame(3,2);
t30 = t50 ^ 2;
t28 = t47 ^ 2;
t26 = t44 ^ 2;
t21 = cos(t33);
t20 = cos(t32);
t19 = cos(t31);
t18 = sin(t33);
t17 = sin(t32);
t16 = sin(t31);
t3 = t81 * t85 + (t29 - 0.1e1) * t42 * pkin(2);
t2 = t82 * t87 + (t27 - 0.1e1) * t39 * pkin(2);
t1 = t83 * t89 + (t25 - 0.1e1) * t36 * pkin(2);
t4 = [((t18 * t90 + t21 * t93) * t30 + (t18 * t84 - t63 * t21) * t50 - t3 * t21 - t18 * t69) * t72, ((-t18 * t93 + t21 * t90) * t30 + (t63 * t18 + t21 * t84) * t50 + t3 * t18 - t21 * t69) * t72, ((sin(-qJ(2,1) + t80) + sin(qJ(2,1) + t79)) * t100 + (-cos(t61 - 0.2e1 * qJ(3,1) + qJ(1,1)) - cos(t60 + t59 + qJ(1,1)) - 0.2e1 * t51) * pkin(2) + (-cos(t61 + t80) - cos(t60 + t79) - cos(t80) - cos(t79)) * pkin(1)) / ((-sin(t59 + qJ(2,1)) + t41) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(qJ(2,1) + qJ(3,1))) * pkin(1)) * t99; ((t17 * t91 + t20 * t94) * t28 + (t17 * t86 - t64 * t20) * t47 - t2 * t20 - t17 * t70) * t73, ((-t17 * t94 + t20 * t91) * t28 + (t64 * t17 + t20 * t86) * t47 + t2 * t17 - t20 * t70) * t73, ((sin(-qJ(2,2) + t78) + sin(qJ(2,2) + t77)) * t100 + (-cos(t58 - 0.2e1 * qJ(3,2) + qJ(1,2)) - cos(t57 + t56 + qJ(1,2)) - 0.2e1 * t48) * pkin(2) + (-cos(t58 + t78) - cos(t57 + t77) - cos(t78) - cos(t77)) * pkin(1)) / ((-sin(t56 + qJ(2,2)) + t38) * pkin(2) + (-sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2))) * pkin(1)) * t99; ((t16 * t92 + t19 * t95) * t26 + (t16 * t88 - t65 * t19) * t44 - t1 * t19 - t16 * t71) * t74, ((-t16 * t95 + t19 * t92) * t26 + (t65 * t16 + t19 * t88) * t44 + t1 * t16 - t19 * t71) * t74, ((sin(-qJ(2,3) + t76) + sin(qJ(2,3) + t75)) * t100 + (-cos(t55 - 0.2e1 * qJ(3,3) + qJ(1,3)) - cos(t54 + t53 + qJ(1,3)) - 0.2e1 * t45) * pkin(2) + (-cos(t55 + t76) - cos(t54 + t75) - cos(t76) - cos(t75)) * pkin(1)) / ((-sin(t53 + qJ(2,3)) + t35) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(qJ(2,3) + qJ(3,3))) * pkin(1)) * t99;];
Jinv  = t4;
