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
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:14
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR12V1G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3P3A1_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:13:52
% EndTime: 2020-08-06 19:13:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->19), mult. (60->45), div. (9->3), fcn. (60->18), ass. (0->35)
t17 = sin(qJ(1,3));
t23 = cos(qJ(1,3));
t16 = sin(qJ(2,3));
t22 = cos(qJ(2,3));
t28 = pkin(1) + pkin(2);
t32 = t16 * qJ(3,3) + t28 * t22;
t40 = -t17 * pkin(4) + t32 * t23;
t19 = sin(qJ(1,2));
t25 = cos(qJ(1,2));
t18 = sin(qJ(2,2));
t24 = cos(qJ(2,2));
t33 = t18 * qJ(3,2) + t28 * t24;
t39 = -t19 * pkin(4) + t33 * t25;
t21 = sin(qJ(1,1));
t27 = cos(qJ(1,1));
t20 = sin(qJ(2,1));
t26 = cos(qJ(2,1));
t34 = t20 * qJ(3,1) + t28 * t26;
t38 = -t21 * pkin(4) + t34 * t27;
t31 = 0.1e1 / qJ(3,1);
t30 = 0.1e1 / qJ(3,2);
t29 = 0.1e1 / qJ(3,3);
t15 = legFrame(1,2);
t14 = legFrame(2,2);
t13 = legFrame(3,2);
t12 = cos(t15);
t11 = cos(t14);
t10 = cos(t13);
t9 = sin(t15);
t8 = sin(t14);
t7 = sin(t13);
t3 = -t26 * qJ(3,1) + t28 * t20;
t2 = -t24 * qJ(3,2) + t28 * t18;
t1 = -t22 * qJ(3,3) + t28 * t16;
t4 = [(t38 * t12 + t3 * t9) * t31, (t12 * t3 - t38 * t9) * t31, (-t27 * pkin(4) - t21 * t34) * t31; (t39 * t11 + t2 * t8) * t30, (t11 * t2 - t39 * t8) * t30, (-t25 * pkin(4) - t19 * t33) * t30; (t1 * t7 + t40 * t10) * t29, (t10 * t1 - t40 * t7) * t29, (-t23 * pkin(4) - t17 * t32) * t29;];
Jinv  = t4;
