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

function Jinv = P3RRRRR1V1G3P3A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G3P3A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G3P3A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G3P3A3_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G3P3A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G3P3A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:38:46
% EndTime: 2020-08-07 03:38:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (36->15), mult. (81->44), div. (27->5), fcn. (90->24), ass. (0->38)
t22 = sin(qJ(3,3));
t45 = pkin(3) * t22;
t24 = sin(qJ(3,2));
t44 = pkin(3) * t24;
t26 = sin(qJ(3,1));
t43 = pkin(3) * t26;
t23 = sin(qJ(2,3));
t28 = cos(qJ(2,3));
t7 = pkin(3) * cos(qJ(3,3)) + pkin(2);
t4 = -t23 * t45 + t7 * t28;
t42 = cos(qJ(1,3)) * t4;
t25 = sin(qJ(2,2));
t30 = cos(qJ(2,2));
t8 = pkin(3) * cos(qJ(3,2)) + pkin(2);
t5 = -t25 * t44 + t8 * t30;
t41 = cos(qJ(1,2)) * t5;
t27 = sin(qJ(2,1));
t32 = cos(qJ(2,1));
t9 = pkin(3) * cos(qJ(3,1)) + pkin(2);
t6 = -t27 * t43 + t9 * t32;
t40 = cos(qJ(1,1)) * t6;
t39 = 0.1e1 / pkin(2) / pkin(3);
t38 = 0.1e1 / t22 * t39;
t37 = 0.1e1 / t24 * t39;
t36 = 0.1e1 / t26 * t39;
t21 = legFrame(1,2);
t20 = legFrame(2,2);
t19 = legFrame(3,2);
t15 = cos(t21);
t14 = cos(t20);
t13 = cos(t19);
t12 = sin(t21);
t11 = sin(t20);
t10 = sin(t19);
t3 = t27 * t9 + t32 * t43;
t2 = t25 * t8 + t30 * t44;
t1 = t23 * t7 + t28 * t45;
t16 = [(t12 * t3 - t15 * t40) * t36, (t12 * t40 + t3 * t15) * t36, sin(qJ(1,1)) * t6 * t36; (t11 * t2 - t14 * t41) * t37, (t11 * t41 + t2 * t14) * t37, sin(qJ(1,2)) * t5 * t37; (t10 * t1 - t13 * t42) * t38, (t1 * t13 + t10 * t42) * t38, sin(qJ(1,3)) * t4 * t38;];
Jinv  = t16;
