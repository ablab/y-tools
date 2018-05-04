#include <gtest/gtest.h>
#include <logger/log_writers.hpp>
#include "../io/Record.hpp"
#include "../ig_tools/utils/string_tools.hpp"
#include "../io/RecordIO.hpp"

const std::string TEST_DATA_PATH = "test_dataset/sequence_io_test_data.fastq";

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

class SequenceIOTest : public ::testing::Test {
public:
    void SetUp() override { create_console_logger(); }
};

TEST_F(SequenceIOTest, CheckRecordCreation) {
    using namespace YTools::IO;
    const auto& int_field = RecordField<size_t>("int_field");
    const auto& string_field = RecordField<std::string>("string_field");

    auto record1 = internal::Record::Create().Set("int_field", static_cast<size_t>(129)).Set("string_field", static_cast<std::string>("string_value")).Build();

    ASSERT_EQ(129, record1.GetField<size_t>("int_field"));
    ASSERT_EQ(129, record1.GetField(int_field));
    ASSERT_EQ("string_value", record1.GetField<std::string>("string_field"));
    ASSERT_EQ("string_value", record1.GetField(string_field));

    auto record2 = internal::Record::Create()
            .Set(SEQUENCE, seqan::Dna5String("ACGT"))
            .Set(CLUSTER_ID, static_cast<std::string>("my_cluster"))
            .Set(CLUSTER_SIZE, 100500ul).Build();

    ASSERT_EQ(3, record2.GetFieldNames().size());
    ASSERT_EQ("ACGT", seqan_string_to_string(record2.GetField<seqan::Dna5String>(SEQUENCE.GetName())));
    ASSERT_EQ("ACGT", seqan_string_to_string(record2.GetField(SEQUENCE)));
    ASSERT_EQ("my_cluster", record2.GetField<std::string>(CLUSTER_ID.GetName()));
    ASSERT_EQ("my_cluster", record2.GetField(CLUSTER_ID));
    ASSERT_EQ(100500ul, record2.GetField<size_t>(CLUSTER_SIZE.GetName()));
    ASSERT_EQ(100500ul, record2.GetField(CLUSTER_SIZE));
}

TEST_F(SequenceIOTest, CheckGetWrongTypeFromRecord) {
    using namespace YTools::IO;
    auto record = internal::Record::Create().Set("int_field", 129).Build();
    ASSERT_ANY_THROW(record.GetField<std::string>("int_field"));
}

TEST_F(SequenceIOTest, CheckGetAbsentFieldFromRecord) {
    using namespace YTools::IO;
    auto record = internal::Record::Create().Set("int_field", 129).Build();
    ASSERT_ANY_THROW(record.GetField<std::string>("string_field"));
}

TEST_F(SequenceIOTest, CheckIgrecOldClusterSizeIO) {
    using namespace YTools::IO;

    path::CheckFileExistenceFATAL(TEST_DATA_PATH);
    std::ifstream test_data_stream(TEST_DATA_PATH);

    IgrecClusterSizeIO seq_io;
    const auto& records = seq_io.ReadRecords(test_data_stream);
    ASSERT_EQ(2, records.size());
    ASSERT_EQ(4, records[0].GetFieldNames().size());
    ASSERT_EQ("0", records[0].GetField(CLUSTER_ID));
    ASSERT_EQ(24, records[0].GetField(CLUSTER_SIZE));
    ASSERT_EQ("ACGTACGTACGTACGTACGTA", seqan_string_to_string(records[0].GetField(SEQUENCE)));
    ASSERT_EQ("ACGT!\"#$%&xyz{|}~ghij", seqan_string_to_string(records[0].GetField(QUALITY)));
    ASSERT_EQ(4, records[1].GetFieldNames().size());
    ASSERT_EQ("my_cluster", records[1].GetField(CLUSTER_ID));
    ASSERT_EQ(3, records[1].GetField(CLUSTER_SIZE));
    ASSERT_EQ("ACGTACGT", seqan_string_to_string(records[1].GetField(SEQUENCE)));
    ASSERT_EQ("!!!!!!!!", seqan_string_to_string(records[1].GetField(QUALITY)));

    std::stringstream serialized_stream;
    seq_io.WriteRecords(serialized_stream, records);
    auto ifs = std::ifstream(TEST_DATA_PATH);
    std::string test_file_content(std::istreambuf_iterator<char>(ifs), {});
    ASSERT_EQ(test_file_content, serialized_stream.str());
}
