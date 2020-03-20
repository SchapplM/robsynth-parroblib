% Return the minimum parameter vector for
% P6PRRRRR6V2G1P4A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-13 01:01
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P6PRRRRR6V2G1P4A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(6,3)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PRRRRR6V2G1P4A0_minimal_parameter_para: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6PRRRRR6V2G1P4A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6PRRRRR6V2G1P4A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P6PRRRRR6V2G1P4A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PRRRRR6V2G1P4A0_minimal_parameter_para: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
unknown=NaN(21,1);
t2 = sin(pkin(5));
t3 = t2 ^ 2;
t8 = pkin(2) ^ 2;
t9 = pkin(8) ^ 2;
unknown(1,1) = m(1) + m(2) + m(3);
unknown(2,1) = Ifges(2,3) + Ifges(3,2) * t3 + 0.2e1 * mrSges(3,3) * t3 * pkin(8) + m(3) * (t3 * t9 + t8);
unknown(3,1) = m(3) * pkin(2) + mrSges(2,1);
unknown(4,1) = -m(3) * t2 * pkin(8) - mrSges(3,3) * t2 + mrSges(2,2);
unknown(5,1) = Ifges(3,1) - Ifges(3,2);
unknown(6,1) = Ifges(3,4);
unknown(7,1) = Ifges(3,5);
unknown(8,1) = Ifges(3,6);
unknown(9,1) = Ifges(3,3);
unknown(10,1) = mrSges(3,1);
unknown(11,1) = mrSges(3,2);
unknown(12,1) = Ifges(4,1);
unknown(13,1) = Ifges(4,4);
unknown(14,1) = Ifges(4,5);
unknown(15,1) = Ifges(4,2);
unknown(16,1) = Ifges(4,6);
unknown(17,1) = Ifges(4,3);
unknown(18,1) = mrSges(4,1);
unknown(19,1) = mrSges(4,2);
unknown(20,1) = mrSges(4,3);
unknown(21,1) = m(4);
MPV  = unknown;
