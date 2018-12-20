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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
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
% Datum: 2018-12-20 17:40
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3PRP2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A1_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:40:11
% EndTime: 2018-12-20 17:40:11
% DurationCPUTime: 0.34s
% Computational Cost: add. (243->88), mult. (492->139), div. (9->3), fcn. (150->14), ass. (0->81)
t61 = (pkin(2) ^ 2);
t86 = t61 + 1;
t85 = (pkin(2) * qJ(3,1));
t84 = (pkin(2) * qJ(3,2));
t83 = (pkin(2) * qJ(3,3));
t48 = cos(qJ(2,3));
t82 = sin(qJ(2,3)) * t48;
t49 = cos(qJ(2,2));
t81 = sin(qJ(2,2)) * t49;
t50 = cos(qJ(2,1));
t80 = sin(qJ(2,1)) * t50;
t58 = koppelP(3,1);
t79 = t86 * t58;
t59 = koppelP(2,1);
t78 = t86 * t59;
t60 = koppelP(1,1);
t77 = t86 * t60;
t52 = (qJ(3,3) ^ 2);
t25 = -t52 + t86;
t53 = (qJ(3,2) ^ 2);
t26 = -t53 + t86;
t54 = (qJ(3,1) ^ 2);
t27 = -t54 + t86;
t57 = koppelP(1,2);
t76 = (t57 * t85);
t75 = t60 * t85;
t56 = koppelP(2,2);
t74 = (t56 * t84);
t73 = t59 * t84;
t55 = koppelP(3,2);
t72 = (t55 * t83);
t71 = t58 * t83;
t44 = legFrame(1,3);
t33 = cos(t44);
t70 = t33 * t85;
t43 = legFrame(2,3);
t32 = cos(t43);
t69 = t32 * t84;
t42 = legFrame(3,3);
t31 = cos(t42);
t68 = t31 * t83;
t28 = sin(t42);
t67 = t28 * t83;
t29 = sin(t43);
t66 = t29 * t84;
t30 = sin(t44);
t65 = t30 * t85;
t13 = -t52 * t58 + 2 * t72 + t79;
t14 = t25 * t55 - 2 * t71;
t51 = xP(3);
t34 = sin(t51);
t35 = cos(t51);
t64 = t34 * t13 + t14 * t35;
t15 = -t53 * t59 + 2 * t74 + t78;
t16 = t26 * t56 - 2 * t73;
t63 = t34 * t15 + t16 * t35;
t17 = -t54 * t60 + 2 * t76 + t77;
t18 = t27 * t57 - 2 * t75;
t62 = t34 * t17 + t18 * t35;
t41 = t50 ^ 2;
t40 = t49 ^ 2;
t39 = t48 ^ 2;
t24 = t61 * t57 + t57 - t75;
t23 = t76 + t77;
t22 = t61 * t56 + t56 - t73;
t21 = t74 + t78;
t20 = t61 * t55 + t55 - t71;
t19 = t72 + t79;
t12 = t27 * t30 + 0.2e1 * t70;
t11 = t26 * t29 + 0.2e1 * t69;
t10 = t25 * t28 + 0.2e1 * t68;
t9 = t33 * t27 - 0.2e1 * t65;
t8 = t32 * t26 - 0.2e1 * t66;
t7 = t31 * t25 - 0.2e1 * t67;
t6 = 0.1e1 / (t27 * t41 + 0.2e1 * t80 * t85 - t54 - t86);
t5 = 0.1e1 / (t26 * t40 + 0.2e1 * t81 * t84 - t53 - t86);
t4 = 0.1e1 / (t25 * t39 + 0.2e1 * t82 * t83 - t52 - t86);
t3 = t17 * t35 - t34 * t18;
t2 = t15 * t35 - t34 * t16;
t1 = t13 * t35 - t34 * t14;
t36 = [(t12 * t80 - t61 * t33 + t9 * t41 - t33 + t65) * t6 (t12 * t41 - t30 * t61 - t9 * t80 - t30 - t70) * t6 ((t30 * t3 - t62 * t33) * t41 - (t3 * t33 + t62 * t30) * t80 + (t34 * t23 + t24 * t35) * t33 - t30 * (t23 * t35 - t34 * t24)) * t6; (t11 * t81 - t61 * t32 + t8 * t40 - t32 + t66) * t5 (t11 * t40 - t29 * t61 - t8 * t81 - t29 - t69) * t5 ((t29 * t2 - t63 * t32) * t40 - (t2 * t32 + t63 * t29) * t81 + (t34 * t21 + t22 * t35) * t32 - t29 * (t21 * t35 - t34 * t22)) * t5; (t10 * t82 - t61 * t31 + t7 * t39 - t31 + t67) * t4 (t10 * t39 - t28 * t61 - t7 * t82 - t28 - t68) * t4 ((t28 * t1 - t64 * t31) * t39 - (t1 * t31 + t64 * t28) * t82 + (t34 * t19 + t20 * t35) * t31 - t28 * (t19 * t35 - t34 * t20)) * t4;];
Jinv  = t36;
