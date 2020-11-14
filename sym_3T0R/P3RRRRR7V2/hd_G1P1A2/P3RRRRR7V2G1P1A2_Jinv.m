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
% Datum: 2020-08-07 09:40
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G1P1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G1P1A2_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:40:14
% EndTime: 2020-08-07 09:40:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (303->106), mult. (222->70), div. (12->7), fcn. (102->69), ass. (0->58)
t78 = 2 * pkin(1);
t46 = (-pkin(7) - pkin(6) - pkin(5));
t77 = -2 * t46;
t75 = 2 * t46;
t40 = legFrame(3,3) + qJ(1,3);
t34 = qJ(3,3) + t40;
t17 = qJ(2,3) + t34;
t35 = -qJ(3,3) + t40;
t18 = -qJ(2,3) + t35;
t74 = sin(t17) + sin(t18);
t41 = legFrame(2,3) + qJ(1,2);
t36 = qJ(3,2) + t41;
t23 = qJ(2,2) + t36;
t37 = -qJ(3,2) + t41;
t24 = -qJ(2,2) + t37;
t73 = sin(t23) + sin(t24);
t42 = legFrame(1,3) + qJ(1,1);
t38 = qJ(3,1) + t42;
t29 = qJ(2,1) + t38;
t39 = -qJ(3,1) + t42;
t30 = -qJ(2,1) + t39;
t72 = sin(t29) + sin(t30);
t43 = sin(qJ(2,3) + qJ(3,3));
t50 = sin(qJ(3,3));
t54 = 0.2e1 * qJ(3,3);
t63 = pkin(2) ^ 2;
t71 = 0.1e1 / (-t63 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (t50 * t78 + pkin(2) * t43 + (-sin(qJ(2,3)) + sin(t54 + qJ(2,3))) * pkin(3))) / 0.2e1;
t44 = sin(qJ(2,2) + qJ(3,2));
t51 = sin(qJ(3,2));
t57 = 0.2e1 * qJ(3,2);
t70 = 0.1e1 / (-t63 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (t51 * t78 + pkin(2) * t44 + (sin(t57 + qJ(2,2)) - sin(qJ(2,2))) * pkin(3))) / 0.2e1;
t45 = sin(qJ(2,1) + qJ(3,1));
t52 = sin(qJ(3,1));
t60 = 0.2e1 * qJ(3,1);
t69 = 0.1e1 / (-t63 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (t52 * t78 + pkin(2) * t45 + (sin(t60 + qJ(2,1)) - sin(qJ(2,1))) * pkin(3))) / 0.2e1;
t68 = cos(t17) + cos(t18);
t67 = cos(t23) + cos(t24);
t66 = cos(t29) + cos(t30);
t64 = 0.1e1 / pkin(2);
t62 = -0.2e1 * qJ(2,1);
t61 = 0.2e1 * qJ(2,1);
t59 = -0.2e1 * qJ(2,2);
t58 = 0.2e1 * qJ(2,2);
t56 = -0.2e1 * qJ(2,3);
t55 = 0.2e1 * qJ(2,3);
t33 = t62 + t39;
t32 = t61 + t38;
t31 = -0.2e1 * qJ(3,1) + t62 + t42;
t28 = t60 + t61 + t42;
t27 = t59 + t37;
t26 = t58 + t36;
t25 = -0.2e1 * qJ(3,2) + t59 + t41;
t22 = t57 + t58 + t41;
t21 = t56 + t35;
t20 = t55 + t34;
t19 = -0.2e1 * qJ(3,3) + t56 + t40;
t16 = t54 + t55 + t40;
t1 = [(t66 * t78 + t72 * t77 + (cos(t31) + cos(t28) + 0.2e1 * cos(t42)) * pkin(3) + (cos(t33) + cos(t32) + cos(t39) + cos(t38)) * pkin(2)) * t69, (t72 * t78 + t66 * t75 + (sin(t31) + sin(t28) + 0.2e1 * sin(t42)) * pkin(3) + (sin(t33) + sin(t32) + sin(t39) + sin(t38)) * pkin(2)) * t69, t64 * t45 / t52; (t67 * t78 + t73 * t77 + (cos(t25) + cos(t22) + 0.2e1 * cos(t41)) * pkin(3) + (cos(t27) + cos(t26) + cos(t37) + cos(t36)) * pkin(2)) * t70, (t73 * t78 + t67 * t75 + (sin(t25) + sin(t22) + 0.2e1 * sin(t41)) * pkin(3) + (sin(t27) + sin(t26) + sin(t37) + sin(t36)) * pkin(2)) * t70, t64 * t44 / t51; (t68 * t78 + t74 * t77 + (cos(t19) + cos(t16) + 0.2e1 * cos(t40)) * pkin(3) + (cos(t21) + cos(t20) + cos(t35) + cos(t34)) * pkin(2)) * t71, (t74 * t78 + t68 * t75 + (sin(t19) + sin(t16) + 0.2e1 * sin(t40)) * pkin(3) + (sin(t21) + sin(t20) + sin(t35) + sin(t34)) * pkin(2)) * t71, t64 * t43 / t50;];
Jinv  = t1;
