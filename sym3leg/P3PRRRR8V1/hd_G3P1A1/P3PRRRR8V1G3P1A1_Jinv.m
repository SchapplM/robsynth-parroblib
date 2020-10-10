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
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:15
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRRRR8V1G3P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3P1A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:15:24
% EndTime: 2020-08-06 17:15:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->24), mult. (210->63), div. (9->3), fcn. (210->22), ass. (0->44)
t49 = pkin(2) * sin(qJ(3,3));
t48 = pkin(2) * sin(qJ(3,2));
t47 = pkin(2) * sin(qJ(3,1));
t46 = pkin(2) * cos(qJ(3,3));
t45 = pkin(2) * cos(qJ(3,2));
t44 = pkin(2) * cos(qJ(3,1));
t30 = sin(qJ(2,3));
t36 = cos(qJ(2,3));
t10 = -t36 * pkin(5) + t30 * t46;
t23 = sin(pkin(3));
t25 = cos(pkin(3));
t43 = -t10 * t25 + t23 * t49;
t32 = sin(qJ(2,2));
t38 = cos(qJ(2,2));
t11 = -t38 * pkin(5) + t32 * t45;
t42 = -t11 * t25 + t23 * t48;
t34 = sin(qJ(2,1));
t40 = cos(qJ(2,1));
t12 = -t40 * pkin(5) + t34 * t44;
t41 = -t12 * t25 + t23 * t47;
t28 = legFrame(1,2);
t27 = legFrame(2,2);
t26 = legFrame(3,2);
t24 = cos(pkin(6));
t22 = sin(pkin(6));
t21 = cos(t28);
t20 = cos(t27);
t19 = cos(t26);
t18 = sin(t28);
t17 = sin(t27);
t16 = sin(t26);
t15 = pkin(5) * t34 + t40 * t44;
t14 = pkin(5) * t32 + t38 * t45;
t13 = pkin(5) * t30 + t36 * t46;
t9 = t12 * t23 + t25 * t47;
t8 = t11 * t23 + t25 * t48;
t7 = t10 * t23 + t25 * t49;
t6 = 0.1e1 / t9;
t5 = 0.1e1 / t8;
t4 = 0.1e1 / t7;
t3 = t24 * t15 + t41 * t22;
t2 = t24 * t14 + t42 * t22;
t1 = t24 * t13 + t43 * t22;
t29 = [(t9 * t18 + t3 * t21) * t6, (-t3 * t18 + t21 * t9) * t6, (-t22 * t15 + t41 * t24) * t6; (t8 * t17 + t2 * t20) * t5, (-t2 * t17 + t20 * t8) * t5, (-t22 * t14 + t42 * t24) * t5; (t1 * t19 + t7 * t16) * t4, (-t1 * t16 + t19 * t7) * t4, (-t22 * t13 + t43 * t24) * t4;];
Jinv  = t29;
