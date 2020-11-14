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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-07 08:51
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V1G2P2A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2P2A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2P2A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2P2A3_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2P2A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2P2A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:49:55
% EndTime: 2020-08-07 08:49:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (54->22), mult. (96->50), div. (27->5), fcn. (105->24), ass. (0->45)
t22 = sin(qJ(3,3));
t49 = pkin(2) * t22;
t25 = sin(qJ(3,2));
t48 = pkin(2) * t25;
t28 = sin(qJ(3,1));
t47 = pkin(2) * t28;
t46 = 0.1e1 / pkin(1) / pkin(2);
t45 = 0.1e1 / t22 * t46;
t44 = 0.1e1 / t25 * t46;
t43 = 0.1e1 / t28 * t46;
t23 = sin(qJ(2,3));
t31 = cos(qJ(2,3));
t7 = cos(qJ(3,3)) * pkin(2) + pkin(1);
t42 = t23 * t49 - t7 * t31;
t26 = sin(qJ(2,2));
t33 = cos(qJ(2,2));
t8 = cos(qJ(3,2)) * pkin(2) + pkin(1);
t41 = t26 * t48 - t8 * t33;
t29 = sin(qJ(2,1));
t35 = cos(qJ(2,1));
t9 = cos(qJ(3,1)) * pkin(2) + pkin(1);
t40 = t29 * t47 - t9 * t35;
t37 = pkin(4) + pkin(5);
t36 = cos(qJ(1,1));
t34 = cos(qJ(1,2));
t32 = cos(qJ(1,3));
t30 = sin(qJ(1,1));
t27 = sin(qJ(1,2));
t24 = sin(qJ(1,3));
t21 = legFrame(1,2);
t20 = legFrame(2,2);
t19 = legFrame(3,2);
t15 = cos(t21);
t14 = cos(t20);
t13 = cos(t19);
t12 = sin(t21);
t11 = sin(t20);
t10 = sin(t19);
t6 = t29 * t9 + t35 * t47;
t5 = t26 * t8 + t33 * t48;
t4 = t23 * t7 + t31 * t49;
t3 = -t40 * t30 - t36 * t37;
t2 = -t41 * t27 - t34 * t37;
t1 = -t42 * t24 - t32 * t37;
t16 = [(-t6 * t12 - t3 * t15) * t43, (t3 * t12 - t15 * t6) * t43, (-t30 * t37 + t40 * t36) * t43; (-t5 * t11 - t2 * t14) * t44, (t2 * t11 - t14 * t5) * t44, (-t27 * t37 + t41 * t34) * t44; (-t1 * t13 - t4 * t10) * t45, (t1 * t10 - t13 * t4) * t45, (-t24 * t37 + t42 * t32) * t45;];
Jinv  = t16;
