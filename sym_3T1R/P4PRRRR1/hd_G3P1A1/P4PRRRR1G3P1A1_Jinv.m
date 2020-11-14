% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [4x4]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:22
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4PRRRR1G3P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(2,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:16:48
% EndTime: 2020-03-02 19:16:48
% DurationCPUTime: 0.20s
% Computational Cost: add. (28->20), mult. (56->56), div. (24->8), fcn. (108->26), ass. (0->44)
t43 = koppelP(1,1);
t42 = koppelP(2,1);
t41 = koppelP(3,1);
t40 = koppelP(4,1);
t39 = koppelP(1,2);
t38 = koppelP(2,2);
t37 = koppelP(3,2);
t36 = koppelP(4,2);
t35 = xP(4);
t34 = cos(qJ(2,1));
t33 = cos(qJ(2,2));
t32 = cos(qJ(2,3));
t31 = sin(qJ(2,1));
t30 = sin(qJ(2,2));
t29 = sin(qJ(2,3));
t28 = legFrame(1,2);
t27 = legFrame(2,2);
t26 = legFrame(3,2);
t25 = legFrame(4,2);
t24 = cos(qJ(2,4));
t23 = sin(qJ(2,4));
t22 = 0.1e1 / t31;
t21 = 0.1e1 / t30;
t20 = 0.1e1 / t29;
t19 = 0.1e1 / t23;
t18 = cos(t35);
t17 = sin(t35);
t16 = cos(t28);
t15 = cos(t27);
t14 = cos(t26);
t13 = cos(t25);
t12 = sin(t28);
t11 = sin(t27);
t10 = sin(t26);
t9 = sin(t25);
t8 = t31 * t12 + t16 * t34;
t7 = -t12 * t34 + t16 * t31;
t6 = t30 * t11 + t15 * t33;
t5 = -t11 * t33 + t15 * t30;
t4 = t29 * t10 + t14 * t32;
t3 = -t10 * t32 + t14 * t29;
t2 = t13 * t24 + t23 * t9;
t1 = t13 * t23 - t9 * t24;
t44 = [t8 * t22, t7 * t22, sin(qJ(3,1)) / cos(qJ(3,1)) * t22, (-t8 * (t17 * t43 + t18 * t39) + t7 * (-t17 * t39 + t18 * t43)) * t22; t6 * t21, t5 * t21, sin(qJ(3,2)) / cos(qJ(3,2)) * t21, (-t6 * (t17 * t42 + t18 * t38) + t5 * (-t17 * t38 + t18 * t42)) * t21; t4 * t20, t3 * t20, sin(qJ(3,3)) / cos(qJ(3,3)) * t20, (-t4 * (t17 * t41 + t18 * t37) + t3 * (-t17 * t37 + t18 * t41)) * t20; t2 * t19, t1 * t19, sin(qJ(3,4)) / cos(qJ(3,4)) * t19, (-t2 * (t17 * t40 + t18 * t36) + t1 * (-t17 * t36 + t18 * t40)) * t19;];
Jinv  = t44;
