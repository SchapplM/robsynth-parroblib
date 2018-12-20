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
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3PRP1A0_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:34:44
% EndTime: 2018-12-20 17:34:44
% DurationCPUTime: 0.33s
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
t55 = koppelP(3,2);
t79 = t86 * t55;
t56 = koppelP(2,2);
t78 = t86 * t56;
t57 = koppelP(1,2);
t77 = t86 * t57;
t52 = (qJ(3,3) ^ 2);
t25 = -t52 + t86;
t53 = (qJ(3,2) ^ 2);
t26 = -t53 + t86;
t54 = (qJ(3,1) ^ 2);
t27 = -t54 + t86;
t76 = t57 * t85;
t60 = koppelP(1,1);
t75 = t60 * t85;
t74 = t56 * t84;
t59 = koppelP(2,1);
t73 = t59 * t84;
t72 = t55 * t83;
t58 = koppelP(3,1);
t71 = t58 * t83;
t44 = legFrame(1,3);
t30 = sin(t44);
t70 = t30 * t85;
t43 = legFrame(2,3);
t29 = sin(t43);
t69 = t29 * t84;
t42 = legFrame(3,3);
t28 = sin(t42);
t68 = t28 * t83;
t31 = cos(t42);
t67 = t31 * t83;
t32 = cos(t43);
t66 = t32 * t84;
t33 = cos(t44);
t65 = t33 * t85;
t64 = -t30 * t27 + 0.2e1 * t65;
t63 = -t29 * t26 + 0.2e1 * t66;
t62 = -t28 * t25 + 0.2e1 * t67;
t51 = xP(3);
t41 = t50 ^ 2;
t40 = t49 ^ 2;
t39 = t48 ^ 2;
t35 = cos(t51);
t34 = sin(t51);
t24 = t75 + t77;
t23 = t61 * t60 + t60 - t76;
t22 = t73 + t78;
t21 = t61 * t59 + t59 - t74;
t20 = t71 + t79;
t19 = t61 * t58 + t58 - t72;
t18 = -t54 * t57 + 2 * t75 + t77;
t17 = t27 * t60 - 2 * t76;
t16 = -t53 * t56 + 2 * t73 + t78;
t15 = t26 * t59 - 2 * t74;
t14 = -t52 * t55 + 2 * t71 + t79;
t13 = t25 * t58 - 2 * t72;
t12 = t27 * t33 + 0.2e1 * t70;
t11 = t26 * t32 + 0.2e1 * t69;
t10 = t25 * t31 + 0.2e1 * t68;
t9 = 0.1e1 / (t27 * t41 + 0.2e1 * t80 * t85 - t54 - t86);
t8 = 0.1e1 / (t26 * t40 + 0.2e1 * t81 * t84 - t53 - t86);
t7 = 0.1e1 / (t25 * t39 + 0.2e1 * t82 * t83 - t52 - t86);
t6 = t34 * t17 + t18 * t35;
t5 = t17 * t35 - t34 * t18;
t4 = t34 * t15 + t16 * t35;
t3 = t15 * t35 - t34 * t16;
t2 = t34 * t13 + t14 * t35;
t1 = t13 * t35 - t34 * t14;
t36 = [(-t12 * t80 + t61 * t30 + t64 * t41 + t30 - t65) * t9 (t12 * t41 - t33 * t61 + t64 * t80 - t33 - t70) * t9 ((t6 * t30 + t5 * t33) * t41 + (-t5 * t30 + t6 * t33) * t80 + (-t23 * t35 + t34 * t24) * t33 - (t34 * t23 + t24 * t35) * t30) * t9; (-t11 * t81 + t61 * t29 + t63 * t40 + t29 - t66) * t8 (t11 * t40 - t32 * t61 + t63 * t81 - t32 - t69) * t8 ((t4 * t29 + t3 * t32) * t40 + (-t3 * t29 + t4 * t32) * t81 + (-t21 * t35 + t34 * t22) * t32 - (t34 * t21 + t22 * t35) * t29) * t8; (-t10 * t82 + t61 * t28 + t62 * t39 + t28 - t67) * t7 (t10 * t39 - t31 * t61 + t62 * t82 - t31 - t68) * t7 ((t1 * t31 + t2 * t28) * t39 + (-t1 * t28 + t2 * t31) * t82 + (-t19 * t35 + t34 * t20) * t31 - (t34 * t19 + t20 * t35) * t28) * t7;];
Jinv  = t36;
