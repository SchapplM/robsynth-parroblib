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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:20
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR12V2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1P1A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:01
% EndTime: 2020-08-06 19:20:01
% DurationCPUTime: 0.17s
% Computational Cost: add. (156->65), mult. (192->114), div. (15->6), fcn. (162->18), ass. (0->63)
t47 = pkin(2) + pkin(3);
t48 = 0.1e1 / qJ(3,3);
t34 = sin(qJ(2,3));
t62 = t34 * qJ(3,3);
t40 = cos(qJ(2,3));
t65 = t47 * t40;
t74 = t48 / (pkin(1) + t62 + t65);
t49 = 0.1e1 / qJ(3,2);
t36 = sin(qJ(2,2));
t61 = t36 * qJ(3,2);
t42 = cos(qJ(2,2));
t64 = t47 * t42;
t73 = t49 / (pkin(1) + t61 + t64);
t50 = 0.1e1 / qJ(3,1);
t38 = sin(qJ(2,1));
t60 = t38 * qJ(3,1);
t44 = cos(qJ(2,1));
t63 = t47 * t44;
t72 = t50 / (pkin(1) + t60 + t63);
t35 = sin(qJ(1,3));
t46 = pkin(5) - pkin(6);
t71 = t35 * t46;
t37 = sin(qJ(1,2));
t70 = t37 * t46;
t39 = sin(qJ(1,1));
t69 = t39 * t46;
t41 = cos(qJ(1,3));
t68 = t46 * t41;
t43 = cos(qJ(1,2));
t67 = t46 * t43;
t45 = cos(qJ(1,1));
t66 = t46 * t45;
t59 = (qJ(3,3) + t47) * (-qJ(3,3) + t47) * t40 ^ 2;
t58 = (qJ(3,2) + t47) * (-qJ(3,2) + t47) * t42 ^ 2;
t57 = (qJ(3,1) + t47) * (-qJ(3,1) + t47) * t44 ^ 2;
t10 = pkin(1) + 0.2e1 * t62;
t56 = t10 * t41 + t71;
t11 = pkin(1) + 0.2e1 * t61;
t55 = t11 * t43 + t70;
t12 = pkin(1) + 0.2e1 * t60;
t54 = t12 * t45 + t69;
t13 = pkin(1) * t34 + qJ(3,3);
t53 = t13 * t41 + t34 * t71;
t14 = pkin(1) * t36 + qJ(3,2);
t52 = t14 * t43 + t36 * t70;
t15 = pkin(1) * t38 + qJ(3,1);
t51 = t15 * t45 + t38 * t69;
t33 = legFrame(1,3);
t32 = legFrame(2,3);
t31 = legFrame(3,3);
t21 = cos(t33);
t20 = cos(t32);
t19 = cos(t31);
t18 = sin(t33);
t17 = sin(t32);
t16 = sin(t31);
t6 = t39 * t12 - t66;
t5 = t37 * t11 - t67;
t4 = t35 * t10 - t68;
t3 = t39 * t15 - t38 * t66;
t2 = t37 * t14 - t36 * t67;
t1 = t35 * t13 - t34 * t68;
t7 = [(-(t18 * t39 - t21 * t45) * t57 - (t18 * t6 - t54 * t21) * t63 - qJ(3,1) * (t3 * t18 - t51 * t21)) * t72, ((t18 * t45 + t21 * t39) * t57 + (t18 * t54 + t6 * t21) * t63 + qJ(3,1) * (t18 * t51 + t3 * t21)) * t72, (-t44 * qJ(3,1) + t47 * t38) * t50; (-(t17 * t37 - t20 * t43) * t58 - (t17 * t5 - t55 * t20) * t64 - qJ(3,2) * (t2 * t17 - t52 * t20)) * t73, ((t17 * t43 + t20 * t37) * t58 + (t17 * t55 + t5 * t20) * t64 + qJ(3,2) * (t17 * t52 + t2 * t20)) * t73, (-t42 * qJ(3,2) + t47 * t36) * t49; (-(t16 * t35 - t19 * t41) * t59 - (t16 * t4 - t56 * t19) * t65 - qJ(3,3) * (t1 * t16 - t53 * t19)) * t74, ((t16 * t41 + t19 * t35) * t59 + (t16 * t56 + t4 * t19) * t65 + qJ(3,3) * (t1 * t19 + t16 * t53)) * t74, (-t40 * qJ(3,3) + t47 * t34) * t48;];
Jinv  = t7;
