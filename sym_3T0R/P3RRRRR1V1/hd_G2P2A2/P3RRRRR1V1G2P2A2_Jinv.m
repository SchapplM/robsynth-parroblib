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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% Datum: 2020-08-07 03:38
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR1V1G2P2A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2P2A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2P2A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2P2A2_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2P2A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2P2A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:38:09
% EndTime: 2020-08-07 03:38:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (33->28), mult. (63->49), div. (18->4), fcn. (99->27), ass. (0->37)
t13 = sin(qJ(3,3));
t14 = sin(qJ(2,3));
t22 = cos(qJ(3,3));
t23 = cos(qJ(2,3));
t41 = sin(qJ(1,3)) * (t13 * t14 - t22 * t23);
t28 = 0.1e1 / pkin(2);
t40 = t28 / t13;
t16 = sin(qJ(3,2));
t39 = t28 / t16;
t19 = sin(qJ(3,1));
t38 = t28 / t19;
t17 = sin(qJ(2,2));
t18 = sin(qJ(1,2));
t37 = t17 * t18;
t25 = cos(qJ(2,2));
t36 = t18 * t25;
t20 = sin(qJ(2,1));
t21 = sin(qJ(1,1));
t35 = t20 * t21;
t27 = cos(qJ(2,1));
t34 = t21 * t27;
t33 = -qJ(2,1) - qJ(3,1);
t32 = -qJ(2,2) - qJ(3,2);
t31 = -qJ(2,3) - qJ(3,3);
t30 = -t13 * t23 - t14 * t22;
t26 = cos(qJ(3,1));
t24 = cos(qJ(3,2));
t12 = legFrame(1,2);
t11 = legFrame(2,2);
t10 = legFrame(3,2);
t6 = cos(t12);
t5 = cos(t11);
t4 = cos(t10);
t3 = sin(t12);
t2 = sin(t11);
t1 = sin(t10);
t7 = [((-t20 * t3 + t6 * t34) * t26 + (-t27 * t3 - t6 * t35) * t19) * t38, ((-t20 * t6 - t3 * t34) * t26 + (-t27 * t6 + t3 * t35) * t19) * t38, (cos(qJ(1,1) + t33) + cos(qJ(1,1) - t33)) * t38 / 0.2e1; ((-t17 * t2 + t5 * t36) * t24 + (-t2 * t25 - t5 * t37) * t16) * t39, ((-t17 * t5 - t2 * t36) * t24 + (t2 * t37 - t25 * t5) * t16) * t39, (cos(qJ(1,2) + t32) + cos(qJ(1,2) - t32)) * t39 / 0.2e1; (t30 * t1 - t4 * t41) * t40, (t1 * t41 + t30 * t4) * t40, (cos(qJ(1,3) + t31) + cos(qJ(1,3) - t31)) * t40 / 0.2e1;];
Jinv  = t7;
