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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:31
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR12V1G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3P3A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:31:09
% EndTime: 2020-08-06 18:31:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (66->37), mult. (93->81), div. (15->6), fcn. (96->18), ass. (0->38)
t20 = sin(qJ(3,3));
t4 = t20 * pkin(3) + qJ(2,3);
t1 = 0.1e1 / t4;
t40 = t1 / t20;
t22 = sin(qJ(3,2));
t5 = t22 * pkin(3) + qJ(2,2);
t2 = 0.1e1 / t5;
t39 = 0.1e1 / t22 * t2;
t24 = sin(qJ(3,1));
t6 = t24 * pkin(3) + qJ(2,1);
t3 = 0.1e1 / t6;
t38 = 0.1e1 / t24 * t3;
t26 = cos(qJ(3,3));
t37 = t26 * qJ(2,3);
t28 = cos(qJ(3,2));
t36 = t28 * qJ(2,2);
t30 = cos(qJ(3,1));
t35 = t30 * qJ(2,1);
t13 = pkin(1) + pkin(5) + pkin(6);
t25 = sin(qJ(1,1));
t31 = cos(qJ(1,1));
t34 = qJ(2,1) * t25 + t13 * t31;
t23 = sin(qJ(1,2));
t29 = cos(qJ(1,2));
t33 = qJ(2,2) * t23 + t13 * t29;
t21 = sin(qJ(1,3));
t27 = cos(qJ(1,3));
t32 = qJ(2,3) * t21 + t13 * t27;
t19 = legFrame(1,2);
t18 = legFrame(2,2);
t17 = legFrame(3,2);
t12 = cos(t19);
t11 = cos(t18);
t10 = cos(t17);
t9 = sin(t19);
t8 = sin(t18);
t7 = sin(t17);
t14 = [(t34 * t12 * t24 + t9 * t35 + (t9 * t30 * t24 + (-t30 ^ 2 + 0.1e1) * t12 * t25) * pkin(3)) * t38, ((t12 * pkin(3) * t30 - t34 * t9) * t24 + t25 * pkin(3) * (t30 - 0.1e1) * (t30 + 0.1e1) * t9 + t12 * t35) * t38, (-t13 * t25 + t6 * t31) * t3; (t33 * t11 * t22 + t8 * t36 + (t8 * t28 * t22 + (-t28 ^ 2 + 0.1e1) * t11 * t23) * pkin(3)) * t39, ((t11 * pkin(3) * t28 - t33 * t8) * t22 + t23 * pkin(3) * (t28 - 0.1e1) * (t28 + 0.1e1) * t8 + t11 * t36) * t39, (-t13 * t23 + t5 * t29) * t2; (t32 * t10 * t20 + t7 * t37 + (t7 * t26 * t20 + (-t26 ^ 2 + 0.1e1) * t10 * t21) * pkin(3)) * t40, ((t10 * pkin(3) * t26 - t32 * t7) * t20 + t21 * pkin(3) * (t26 - 0.1e1) * (t26 + 0.1e1) * t7 + t10 * t37) * t40, (-t13 * t21 + t4 * t27) * t1;];
Jinv  = t14;
